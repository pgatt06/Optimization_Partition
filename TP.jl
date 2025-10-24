# ================CONFIGURATION GUROBI================#
# Déclare le chemin local vers Gurobi et sa licence.
ENV["GUROBI_HOME"] = "/Library/gurobi1203/macos_universal2/"
ENV["GRB_LICENSE_FILE"] = "/Users/p.gatt/Documents/3A/Optimisation_discrete/gurobi.lic"

#================IMPORTATION DES PACKAGES================#
# Installe les dépendances si besoin (ne fait rien si elles sont déjà présentes).
import Pkg
Pkg.add(["JuMP", "Gurobi", "MathOptInterface"])
Pkg.add(["Graphs", "GraphRecipes", "Plots", "Colors", "PrettyTables"])

using Graphs
using GraphRecipes
using Plots
using Colors
using PrettyTables
using JuMP
using Gurobi
const INF = 1.0e12

# ================PARAMETRES================#
# Paramètres généraux utilisés par tous les modèles.
const time_limit = 120.0        # en secondes
const gap_limit  = 1e-4          # 0.01 %
const max_rounds = 50
const tol        = 1e-6          # tolérance pour les coupes

#================FONCTION LECTURE DATA ================#
"""
readWeightedGraph_paper(file) -> (E, W)

Lit une instance au format LOCO-UNICAMP.
Retourne :
- E : matrice d’adjacence (Int, n×n) non orientée (symétrisée)
- W : vecteur poids sommets (Int, n)
"""
function readWeightedGraph_paper(file::String)
    open(file, "r") do io
        data = readlines(io)
        line = split(strip(data[2]))
        n = parse(Int, line[1])
        m = parse(Int, line[2])
        E = zeros(Int, n, n)
        W = zeros(Int, n)
        # Poids
        for i in 1:n
            li = split(strip(data[3 + i]))
            W[i] = parse(Int, li[4])
        end
        # Arêtes (symétrisées)
        for i in 1:m
            le = split(strip(data[4 + n + i]))
            u = parse(Int, le[1]) + 1
            v = parse(Int, le[2]) + 1
            E[u, v] = 1
            E[v, u] = 1
        end
        return E, W
    end
end

# ================METHODE 1 : FLOW (Q1)================#
"""
solve_bcpk_flow(E, W, k) -> (obj, classes)

Formulation à flots avec k sources fictives. Renvoie l'objectif et la partition.
"""
function solve_bcpk_flow(E::Array{Int,2}, W::Vector{Int}, k::Int)
    n  = length(W)
    wG = sum(W)

    # k sources fictives
    sources = collect(n+1:n+k)

    # Arcs (orientés) le long des arêtes non orientées + arcs source->sommet
    arcs = Tuple{Int,Int}[]
    for u in 1:n, v in 1:n
        if u != v && E[u,v] == 1
            push!(arcs, (u, v))
        end
    end
    for si in sources, v in 1:n
        push!(arcs, (si, v))
    end

    A = length(arcs)
    arcs_out = [Int[] for _ in 1:(n+k)]
    arcs_in  = [Int[] for _ in 1:(n+k)]
    for a in 1:A
        (u,v) = arcs[a]
        push!(arcs_out[u], a)
        push!(arcs_in[v], a)
    end

    model = Model(() -> Gurobi.Optimizer())
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", time_limit)
    set_optimizer_attribute(model, "MIPGap",    gap_limit)

    @variable(model, f[1:A] >= 0.0)
    @variable(model, y[1:A], Bin)

    s1 = sources[1]
    @objective(model, Max, sum(f[a] for a in arcs_out[s1]))

    # ordre non décroissant des flux par source
    for i in 1:k-1
        si, sip1 = sources[i], sources[i+1]
        @constraint(model, sum(f[a] for a in arcs_out[si]) <= sum(f[a] for a in arcs_out[sip1]))
    end

    # conservation + consommation W[v] aux sommets réels
    for v in 1:n
        @constraint(model, sum(f[a] for a in arcs_in[v]) - sum(f[a] for a in arcs_out[v]) == W[v])
    end

    # activation par y
    for a in 1:A
        @constraint(model, f[a] <= wG * y[a])
    end

    # une sortie par source + au plus un père par sommet réel
    for si in sources
        @constraint(model, sum(y[a] for a in arcs_out[si]) <= 1)
    end
    for v in 1:n
        @constraint(model, sum(y[a] for a in arcs_in[v]) <= 1)
    end

    optimize!(model)
    if !has_values(model)
        return (NaN, [Int[] for _ in 1:k])
    end

    obj  = objective_value(model)
    yval = value.(y)

    # reconstruire la forêt et les classes
    parent = fill(0, n)
    for v in 1:n
        for a in arcs_in[v]
            if yval[a] > 0.5
                parent[v] = arcs[a][1]
                break
            end
        end
    end
    children = [Int[] for _ in 1:(n+k)]
    for v in 1:n
        p = parent[v]
        p != 0 && push!(children[p], v)
    end
    classes = [Int[] for _ in 1:k]
    for i in 1:k
        si = sources[i]
        stack = copy(children[si])
        while !isempty(stack)
            v = pop!(stack)
            push!(classes[i], v)
            append!(stack, children[v])
        end
        sort!(classes[i])
    end
    return obj, classes
end

# ================OUTIL : MIN-CUT (Dinic)================#
"""
mincut_two_layer(nb_nodes, origins, dests, caps, s, t) -> (flow, side)

Dinic sur graphe résiduel construit à partir d'arêtes (origins -> dests) avec capacités caps.
Retourne le flot max et l'ensemble des sommets côté source atteignables (side).
"""
function mincut_two_layer(nb_nodes::Int,
                          origins::Vector{Int}, dests::Vector{Int}, caps::Vector{Float64},
                          s::Int, t::Int)
    adj = [Int[] for _ in 1:nb_nodes]
    to  = Int[]; rev = Int[]; cap = Float64[]

    function add_edge(u::Int, v::Int, c::Float64)
        push!(to, v);   push!(rev, length(to) + 1); push!(cap, c);   push!(adj[u], length(to))
        push!(to, u);   push!(rev, length(to) - 1); push!(cap, 0.0); push!(adj[v], length(to))
    end

    for e in eachindex(origins)
        add_edge(origins[e], dests[e], caps[e])
    end

    level = fill(-1, nb_nodes)
    itcur = fill(1, nb_nodes)

    function bfs!()
        fill!(level, -1)
        level[s] = 0
        q = Int[s]
        head = 1
        while head <= length(q)
            v = q[head]; head += 1
            for idx in adj[v]
                u = to[idx]
                if cap[idx] > 1e-12 && level[u] < 0
                    level[u] = level[v] + 1
                    push!(q, u)
                end
            end
        end
        return level[t] >= 0
    end

    function dfs!(v::Int, fmax::Float64)
        v == t && return fmax
        for k in itcur[v]:length(adj[v])
            itcur[v] = k
            idx = adj[v][k]
            u   = to[idx]
            if cap[idx] > 1e-12 && level[v] < level[u]
                δ = dfs!(u, min(fmax, cap[idx]))
                if δ > 1e-12
                    cap[idx]      -= δ
                    cap[rev[idx]] += δ
                    return δ
                end
            end
        end
        return 0.0
    end

    flow = 0.0
    while bfs!()
        fill!(itcur, 1)
        while true
            δ = dfs!(s, INF)
            δ < 1e-12 && break
            flow += δ
        end
    end

    seen = falses(nb_nodes)
    stack = [s]; seen[s] = true
    while !isempty(stack)
        v = pop!(stack)
        for idx in adj[v]
            u = to[idx]
            if cap[idx] > 1e-12 && !seen[u]
                seen[u] = true
                push!(stack, u)
            end
        end
    end
    side = findall(identity, seen)
    return flow, side
end

# ================ENUMERATION 4-CYCLES================#
"""
enumerate_4cycles(E) -> Vector{Vector{Int}}

Détecte des 4-cycles simples (a-b-c-d-a) sans chordes dans le graphe non orienté.
"""
function enumerate_4cycles(E::Array{Int,2})
    n = size(E,1)
    adj = [Int[] for _ in 1:n]
    for u in 1:n, v in u+1:n
        if E[u,v] == 1
            push!(adj[u], v); push!(adj[v], u)
        end
    end
    cycles = Set{NTuple{4,Int}}()
    for u in 1:n
        for v in adj[u]
            v <= u && continue
            for w in adj[v]
                w == u && continue
                for x in adj[w]
                    (x == v || x == u) && continue
                    (x in adj[u]) || continue        # ferme le cycle
                    (w in adj[u]) && continue        # pas de corde (u-w)
                    (v in adj[x]) && continue        # pas de corde (v-x)
                    cyc = (u,v,w,x)
                    candidates = NTuple{4,Int}[
                        cyc,(v,w,x,u),(w,x,u,v),(x,u,v,w),
                        (u,x,w,v),(x,w,v,u),(w,v,u,x),(v,u,x,w)
                    ]
                    push!(cycles, minimum(candidates))
                end
            end
        end
    end
    return [collect(c) for c in cycles]
end

# ================CROSS INEQUALITIES================#
"""
add_cross_inequalities!(model, x, x_sol, four_cycles, k; tol) -> Int

Ajoute les inégalités croisées violées :
x[a,i1] + x[c,i1] + x[b,i2] + x[d,i2] ≤ 3, i1≠i2
pour chaque 4-cycle (a-b-c-d).
"""
function add_cross_inequalities!(model::JuMP.Model, x, x_sol::Matrix{Float64},
                                 four_cycles::Vector{Vector{Int}},
                                 k::Int; tol::Float64=1e-9)
    cuts = 0
    for cyc in four_cycles
        a,b,c,d = cyc
        for i1 in 1:k, i2 in 1:k
            i1 == i2 && continue
            lhs_val = x_sol[a,i1] + x_sol[c,i1] + x_sol[b,i2] + x_sol[d,i2]
            if lhs_val > 3.0 + tol
                @constraint(model, x[a,i1] + x[c,i1] + x[b,i2] + x[d,i2] <= 3.0)
                cuts += 1
            end
        end
    end
    return cuts
end

# ================METHODE 2 : CUTS (Q2, sans boost)================#
"""
solve_bcpk_coupes_separation(E, W, k) -> (obj, classes)

Formulation cut-based avec séparation de connectivité (sans boost) :
déclenche si x[u,i] + x[v,i] > 1 + tol.
"""
function solve_bcpk_coupes_separation(E::Array{Int,2}, W::Vector{Int}, k::Int)
    n = length(W)
    start_time = time()

    # couples non-arêtes (u < v)
    non_edges = Tuple{Int,Int}[]
    for u in 1:n, v in u+1:n
        if E[u,v] == 0
            push!(non_edges, (u,v))
        end
    end

    model = Model(() -> Gurobi.Optimizer())
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "MIPGap",    gap_limit)
    set_optimizer_attribute(model, "TimeLimit", time_limit)

    @variable(model, x[1:n, 1:k], Bin)

    # ordre des classes par poids total
    @constraint(model, [i in 1:k-1],
        sum(W[v]*x[v,i] for v in 1:n) <= sum(W[v]*x[v,i+1] for v in 1:n)
    )
    # appartenance unique
    @constraint(model, [v in 1:n], sum(x[v,i] for i in 1:k) <= 1)

    @objective(model, Max, sum(W[v]*x[v,1] for v in 1:n))

    # arêtes non orientées
    edges_undirected = Tuple{Int,Int}[]
    for u in 1:n, v in u+1:n
        if E[u,v] == 1
            push!(edges_undirected, (u,v))
        end
    end

    # pré-construction statique des arcs INF
    base_orig = Int[]; base_dest = Int[]; base_cap = Float64[]
    add_arc!(u,v,c) = (push!(base_orig,u); push!(base_dest,v); push!(base_cap,c))
    for (u,v) in edges_undirected
        add_arc!(u + n, v, INF)
        add_arc!(v + n, u, INF)
    end

    function add_violated_cuts!(x_sol::Matrix{Float64})
        cuts_added = 0
        for i in 1:k
            origins = copy(base_orig); dests = copy(base_dest); caps = copy(base_cap)
            for v in 1:n
                push!(origins, v); push!(dests, v+n); push!(caps, float(x_sol[v,i]))
            end
            for (u,v) in non_edges
                # SANS BOOST : condition sur la somme
                if x_sol[u,i] + x_sol[v,i] > 1.0 + tol
                    s = u; t = v + n
                    value, side = mincut_two_layer(2n, origins, dests, caps, s, t)
                    threshold = x_sol[u,i] + x_sol[v,i] - 1.0
                    if value + 1e-9 < threshold
                        in_side = falses(2n)
                        for z in side; in_side[z] = true; end
                        S = Int[]
                        for z in 1:n
                            if in_side[z] && !in_side[z+n]
                                push!(S, z)
                            end
                        end
                        @constraint(model, x[u,i] + x[v,i] - sum(x[z,i] for z in S) <= 1.0)
                        cuts_added += 1
                    end
                end
            end
        end
        return cuts_added
    end

    round = 0
    while (round < max_rounds) && ((time() - start_time) < time_limit)
        optimize!(model)
        has_values(model) || break
        x_sol = value.(x)
        added = add_violated_cuts!(x_sol)
        round += 1
        added == 0 && break
    end

    optimize!(model)
    if !has_values(model)
        return (NaN, [Int[] for _ in 1:k])
    end

    obj  = objective_value(model)
    xval = value.(x)
    classes = [Int[] for _ in 1:k]
    for i in 1:k
        for v in 1:n
            xval[v,i] > 0.5 && push!(classes[i], v)
        end
        sort!(classes[i])
    end
    return obj, classes
end

# ================METHODE 3 : CUTS (Q3, avec boost + cross)================#
"""
solve_bcpk_coupes_separation_boost(E, W, k) -> (obj, classes)

Identique à la méthode 2, mais **avec boost** (x[u,i] > 0.5 et x[v,i] > 0.5)
et ajout des **cross inequalities** sur 4-cycles.
"""
function solve_bcpk_coupes_separation_boost(E::Array{Int,2}, W::Vector{Int}, k::Int)
    n = length(W)
    start_time = time()

    non_edges = Tuple{Int,Int}[]
    for u in 1:n, v in u+1:n
        E[u,v] == 0 && push!(non_edges, (u,v))
    end

    model = Model(() -> Gurobi.Optimizer())
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "MIPGap",    gap_limit)
    set_optimizer_attribute(model, "TimeLimit", time_limit)

    @variable(model, x[1:n, 1:k], Bin)
    @constraint(model, [i in 1:k-1],
        sum(W[v]*x[v,i] for v in 1:n) <= sum(W[v]*x[v,i+1] for v in 1:n)
    )
    @constraint(model, [v in 1:n], sum(x[v,i] for i in 1:k) <= 1)
    @objective(model, Max, sum(W[v]*x[v,1] for v in 1:n))

    edges_undirected = Tuple{Int,Int}[]
    for u in 1:n, v in u+1:n
        E[u,v] == 1 && push!(edges_undirected, (u,v))
    end

    # 4-cycles pour cross inequalities (pré-calcul)
    four_cycles = enumerate_4cycles(E)

    # pré-construction statique des arcs INF
    base_orig = Int[]; base_dest = Int[]; base_cap = Float64[]
    add_arc!(u,v,c) = (push!(base_orig,u); push!(base_dest,v); push!(base_cap,c))
    for (u,v) in edges_undirected
        add_arc!(u + n, v, INF)
        add_arc!(v + n, u, INF)
    end

    function add_violated_cuts!(x_sol::Matrix{Float64})
        cuts_added = 0

        # 1) Connectivité (avec BOOST)
        for i in 1:k
            origins = copy(base_orig); dests = copy(base_dest); caps = copy(base_cap)
            for v in 1:n
                push!(origins, v); push!(dests, v+n); push!(caps, float(x_sol[v,i]))
            end
            for (u,v) in non_edges
                if (x_sol[u,i] > 0.5 + tol) && (x_sol[v,i] > 0.5 + tol)
                    s = u; t = v + n
                    value, side = mincut_two_layer(2n, origins, dests, caps, s, t)
                    threshold = x_sol[u,i] + x_sol[v,i] - 1.0
                    if value + 1e-9 < threshold
                        in_side = falses(2n)
                        for z in side; in_side[z] = true; end
                        S = Int[]
                        for z in 1:n
                            if in_side[z] && !in_side[z+n]
                                push!(S, z)
                            end
                        end
                        @constraint(model, x[u,i] + x[v,i] - sum(x[z,i] for z in S) <= 1.0)
                        cuts_added += 1
                    end
                end
            end
        end

        # 2) Cross inequalities sur 4-cycles
        cuts_added += add_cross_inequalities!(model, x, x_sol, four_cycles, k; tol=tol)

        return cuts_added
    end

    round = 0
    while (round < max_rounds) && ((time() - start_time) < time_limit)
        optimize!(model)
        has_values(model) || break
        x_sol = value.(x)
        added = add_violated_cuts!(x_sol)
        round += 1
        added == 0 && break
    end

    optimize!(model)
    if !has_values(model)
        return (NaN, [Int[] for _ in 1:k])
    end

    obj  = objective_value(model)
    xval = value.(x)
    classes = [Int[] for _ in 1:k]
    for i in 1:k
        for v in 1:n
            xval[v,i] > 0.5 && push!(classes[i], v)
        end
        sort!(classes[i])
    end
    return obj, classes
end

# =================== FONCTION AFFICHAGE ==================#
function print_solution(title::String, obj::Float64, classes::Vector{Vector{Int}}, W::Vector{Int})
    println(title)
    println("Poids minimal maximisé = ", obj)
    for (i,cl) in enumerate(classes)
        w = sum(W[v] for v in cl)
        println("  Classe $(i) (poids=$w) : ", cl)
    end
end

function plot_partition(E, classes; node_colors = [:red, :blue, :green, :orange, :purple])
    n = size(E,1)
    g = SimpleGraph(n)
    for u in 1:n, v in u+1:n
        E[u,v]==1 && add_edge!(g,u,v)
    end
    color = fill(:gray, n)
    for (i,cl) in enumerate(classes)
        for v in cl
            color[v] = node_colors[mod1(i,length(node_colors))]
        end
    end
    graphplot(g, nodecolor=color, nodesize=0.15, markerstrokewidth=0.5)
end


# ================= PETITS TESTS =================#

function check_connectivity(E, classes)
    n = size(E, 1)
    g = SimpleGraph(n)
    for u in 1:n, v in u+1:n
        E[u,v] == 1 && add_edge!(g, u, v)
    end
    for (i,cl) in enumerate(classes)
        sub, _ = induced_subgraph(g, cl)
        if !is_connected(sub)
            println("Classe $i n’est pas connexe : ", cl)
        else
            println("Classe $i est connexe (|cl|=", length(cl), ")")
        end
    end
end


# =================LANCEMENT DES TESTS ===============#
# Lit l'instance de test puis lance les trois formulations.
# Exemple : adapter le chemin à ton installation

const FNAME = "/Users/p.gatt/Documents/3A/Optimisation_discrete/random/rnd_n20/m30/a/rnd_20_30_a_1.in"
E, W = readWeightedGraph_paper(FNAME)
k = 2

println("=== Méthode 1 : FLOW (Q1) ===")
objF, classesF = solve_bcpk_flow(E, W, k)
print_solution("FLOW", objF, classesF, W)

println("\n=== Méthode 2 : CUTS + séparation (Q2) ===")
objC, classesC = solve_bcpk_coupes_separation(E, W, k)
print_solution("CUTS", objC, classesC, W)

println("\n=== Méthode 3 : CUTS + séparation BOOSTÉE (Q3) ===")
objCB, classesCB = solve_bcpk_coupes_separation_boost(E, W, k)
print_solution("CUTS BOOSTÉE", objCB, classesCB, W)


# PLOTS ET TESTS 
check_connectivity(E, classesF)
plot_partition(E, classesF)