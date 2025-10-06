# ==========================
# Initialisation Gurobi
# ==========================
ENV["GUROBI_HOME"] = "/Library/gurobi1203/macos_universal2/"
ENV["GRB_LICENSE_FILE"] = "/Users/p.gatt/Documents/3A/Optimisation_discrete/gurobi.lic"

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

############# Lecture de l’instance papier #############
function readWeightedGraph_paper(file::String)
    open(file, "r") do io
        data = readlines(io)

        # ligne 2 : n, m
        line = split(strip(data[2]))
        n = parse(Int, line[1])
        m = parse(Int, line[2])

        E = zeros(Int, n, n)   # matrice d’adjacence orientée (0/1)
        W = zeros(Int, n)      # poids sommets

        # lignes 4..(3+n) : poids (4e champ)
        for i in 1:n
            line = split(strip(data[3 + i]))
            W[i] = parse(Int, line[4])
        end

        # lignes (5+n)..(4+n+m) : arcs 0-based -> +1
        for i in 1:m
            line = split(strip(data[4 + n + i]))
            orig = parse(Int, line[1]) + 1
            dest = parse(Int, line[2]) + 1
            E[orig, dest] = 1
        end
        return E, W
    end
end

############# 1) Formulation FLOW (Q1) #############
function solve_bcpk_flow(E::Array{Int,2}, W::Vector{Int}, k::Int)
    n  = length(W)
    wG = sum(W)
    sources = collect(n+1:n+k)

    # arcs autorisés : on “symétrise” E (si u-v ou v-u, on autorise u->v)
    arcs = Tuple{Int,Int}[]
    for u in 1:n, v in 1:n
        if u != v && (E[u,v] == 1 || E[v,u] == 1)
            push!(arcs, (u, v))
        end
    end
    # arcs depuis chaque source fictive vers tous les sommets
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

    @variable(model, f[1:A] >= 0.0)
    @variable(model, y[1:A], Bin)

    s1 = sources[1]
    @objective(model, Max, sum(f[a] for a in arcs_out[s1]))

    # ordre croissant des flots sources (≈ équilibre des classes par poids)
    for i in 1:k-1
        si, sip1 = sources[i], sources[i+1]
        @constraint(model, sum(f[a] for a in arcs_out[si]) <= sum(f[a] for a in arcs_out[sip1]))
    end

    # conservation + demande = W[v] sur les sommets réels
    for v in 1:n
        @constraint(model, sum(f[a] for a in arcs_in[v]) - sum(f[a] for a in arcs_out[v]) == W[v])
    end

    # capacité-liens f <= M*y
    for a in 1:A
        @constraint(model, f[a] <= wG * y[a])
    end

    # chaque source connecte au plus 1 arc sortant
    for si in sources
        @constraint(model, sum(y[a] for a in arcs_out[si]) <= 1)
    end
    # chaque sommet a au plus 1 arc entrant sélectionné
    for v in 1:n
        @constraint(model, sum(y[a] for a in arcs_in[v]) <= 1)
    end

    optimize!(model)

    obj  = objective_value(model)
    yval = value.(y)

    # reconstruction parent/enfants/classes
    parent = fill(0, n)
    for v in 1:n
        for a in arcs_in[v]
            if yval[a] > 0.5
                (u, _) = arcs[a]
                parent[v] = u
                break
            end
        end
    end

    children = [Int[] for _ in 1:(n+k)]
    for v in 1:n
        p = parent[v]
        if p != 0
            push!(children[p], v)
        end
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

############# Dinic min-cut sur graphe 2-couches (pour les coupes) #############
# Graphe : 1..n (couche 1), n+1..2n (couche 2)
# Arcs "u2->v1" inf quand (u,v) adjacent ; arcs "v1->v2" capacité = x[v] (LP value ou 0/1).
function mincut_two_layer(nb_nodes::Int,
                          origins::Vector{Int}, dests::Vector{Int}, caps::Vector{Float64},
                          s::Int, t::Int)
    adj = [Int[] for _ in 1:nb_nodes]
    to  = Int[]; rev = Int[]; cap = Float64[]

    function add_edge(u, v, c)
        push!(to, v); push!(rev, length(to)+1); push!(cap, c); push!(adj[u], length(to))
        push!(to, u); push!(rev, length(to)-1); push!(cap, 0.0); push!(adj[v], length(to))
    end
    for e in eachindex(origins)
        add_edge(origins[e], dests[e], caps[e])
    end

    level  = fill(-1, nb_nodes)
    itcur  = fill(1, nb_nodes)

    function bfs!()
        fill!(level, -1); level[s] = 0
        q = Int[s]; head = 1
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

    function dfs!(v, fmax)
        v == t && return fmax
        for k in itcur[v]:length(adj[v])
            itcur[v] = k
            idx = adj[v][k]
            u = to[idx]
            if cap[idx] > 1e-12 && level[v] < level[u]
                δ = dfs!(u, min(fmax, cap[idx]))
                if δ > 1e-12
                    cap[idx] -= δ
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
            pushed = dfs!(s, INF)
            pushed <= 1e-12 && break
            flow += pushed
        end
    end

    # composante atteignable depuis s dans le résiduel
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

############# 2) COUPES – boucle externe (zéro callback) #############
function solve_bcpk_coupes_no_callback(E::Array{Int,2}, W::Vector{Int}, k::Int;
                                       tol::Float64 = 1e-6, max_rounds::Int = 50)
    n = length(W)
    @assert size(E,1) == n && size(E,2) == n

    # Paires non-arêtes
    non_edges = Tuple{Int,Int}[]
    for u in 1:n, v in u+1:n
        if E[u,v] == 0 && E[v,u] == 0
            push!(non_edges, (u,v))
        end
    end

    # Modèle de base
    model = Model(() -> Gurobi.Optimizer())
    set_optimizer_attribute(model, "OutputFlag", 0)

    @variable(model, x[1:n, 1:k], Bin)

    # équilibre de poids croissant entre classes
    @constraint(model, [i in 1:k-1],
        sum(W[v]*x[v,i] for v in 1:n) <= sum(W[v]*x[v,i+1] for v in 1:n)
    )
    # affectation unique
    @constraint(model, [v in 1:n], sum(x[v,i] for i in 1:k) <= 1)

    # objectif : max poids classe 1
    @objective(model, Max, sum(W[v]*x[v,1] for v in 1:n))

    # Prépare arêtes pour le graphe 2-couches
    edges_undirected = Tuple{Int,Int}[]
    for u in 1:n, v in u+1:n
        if E[u,v] == 1 || E[v,u] == 1
            push!(edges_undirected, (u,v))
        end
    end

    # fonction qui construit et ajoute les coupes violées à partir d'une solution binaire
    function add_violated_cuts!(x_sol::Matrix{Float64})
        cuts_added = 0
        for i in 1:k
            # construire la liste d'arcs du graphe 2-couches
            origins = Int[]; dests = Int[]; caps = Float64[]
            add_arc(u::Int, v::Int, c::Float64) = (push!(origins,u); push!(dests,v); push!(caps,c))

            # Arêtes du graphe: u2->v1 et v2->u1 avec capacité INF
            for (u,v) in edges_undirected
                add_arc(u + n, v, INF)
                add_arc(v + n, u, INF)
            end
            # arcs (v1->v2) capacité = x[v,i]
            for v in 1:n
                add_arc(v, v + n, float(x_sol[v,i]))
            end

            # Chercher violations : u,v non-arête et tous deux sélectionnés
            for (u,v) in non_edges
                if x_sol[u,i] > 0.5 + tol && x_sol[v,i] > 0.5 + tol
                    s = u       # u1
                    t = v + n   # v2
                    value, side = mincut_two_layer(2n, origins, dests, caps, s, t)
                    # seuil = x[u,i] + x[v,i] - 1 ; si min-cut < seuil -> coupe violée
                    threshold = x_sol[u,i] + x_sol[v,i] - 1.0
                    if value + 1e-9 < threshold
                        in_side = falses(2n)
                        for z in side
                            in_side[z] = true
                        end
                        S = Int[]
                        for z in 1:n
                            if in_side[z] && !in_side[z+n]
                                push!(S, z)
                            end
                        end
                        # coupe : x_u + x_v - sum_{z∈S} x_z ≤ 1
                        @constraint(model, x[u,i] + x[v,i] - sum(x[z,i] for z in S) <= 1.0)
                        cuts_added += 1
                    end
                end
            end
        end
        return cuts_added
    end

    # Boucle cutting-plane
    for round in 1:max_rounds
        optimize!(model)
        term = termination_status(model)
        if term != MOI.OPTIMAL && term != MOI.TIME_LIMIT && term != MOI.MEMORY_LIMIT
            break
        end

        x_sol = value.(x)
        added = add_violated_cuts!(x_sol)
        # println("Round $round — coupes ajoutées: $added")
        added == 0 && break
    end

    obj = objective_value(model)
    xval = value.(x)

    classes = [Int[] for _ in 1:k]
    for i in 1:k
        for v in 1:n
            if xval[v,i] > 0.5
                push!(classes[i], v)
            end
        end
        sort!(classes[i])
    end
    return obj, classes
end

############# Démo minimale #############
const FNAME = "instance_bcpk.txt"
instance_text = """
nnodes nedges type
8 10 graph
nodename posx posy weight
0 0 0 2
1 0 0 1
2 0 0 1
3 0 0 3
4 0 0 1
5 0 0 3
6 0 0 4
7 0 0 2
endpoint1 endpoint2
0 1
1 2
0 2
2 3
3 4
3 6
4 6
4 5
5 7
6 7
"""
open(FNAME, "w") do io
    write(io, instance_text)
end

E, W = readWeightedGraph_paper(FNAME)
k = 2

println("=== Méthode 1 : FLOW (Q1) ===")
objF, classesF = solve_bcpk_flow(E, W, k)
println("Poids minimal maximisé = ", objF)
for i in 1:k
    println("Classe $i : ", classesF[i])
end

println("\n=== Méthode 2 : COUPES (sans callback) ===")
objC, classesC = solve_bcpk_coupes_no_callback(E, W, k; max_rounds=50)
println("Poids minimal maximisé = ", objC)
for i in 1:k
    println("Classe $i : ", classesC[i])
end
