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
import MathOptInterface
const MOI = MathOptInterface

# ==========================
# Lecture de l’instance
# ==========================
function readWeightedGraph_paper(file::String)

    #Open file, extract data, close file
    
    o = open(file,"r")
    
    data = readlines(o)
    
    close(o)
    
    
    #nombre de sommets et d'aretes
    
    line = split(data[2]," ") #lecture de la ligne numero 2
    
    n=parse(Int64,line[1])
    
    m=parse(Int64,line[2])
    
    
    E = Matrix{Int64}(zeros(n,n))
    
    W = Vector{Int64}(zeros(n))
    
    for i in 1:n
    
    line = split(data[3+i]," ") #lecture de la ligne numero 3+i
    
    poids= parse(Int64,line[4]) # 4eme nombre de la ligne, les autres ne sont pas significatifs
    
    W[i]=poids
    
    end
    
    for i in 1:m
    
    line = split(data[4+n+i]," ") #lecture de la ligne numero 4+n+i
    
    orig = parse(Int64,line[1])
    
    dest = parse(Int64,line[2])
    
    E[orig+1,dest+1] = 1
    
    end
    
    return E, W
    
    end

# ==========================
# 1) Formulation FLOW (Section 4)
# ==========================
function solve_bcpk_flow(E::Array{Int,2}, W::Vector{Int}, k::Int)
    n  = length(W)
    wG = sum(W)
    sources = collect(n+1:n+k)

    arcs = Tuple{Int,Int}[]
    # arcs entre sommets (utilise E symétrisé)
    for u in 1:n, v in 1:n
        if E[u,v] == 1
            push!(arcs, (u,v))
        end
    end
    # arcs source -> sommets
    for si in sources, v in 1:n
        push!(arcs, (si,v))
    end

    A = length(arcs)
    out_arcs = [Int[] for _ in 1:(n+k)]
    in_arcs  = [Int[] for _ in 1:(n+k)]
    for a in 1:A
        (u,v) = arcs[a]
        push!(out_arcs[u], a)
        push!(in_arcs[v], a)
    end

    model = Model(Gurobi.Optimizer)
    set_silent(model)

    @variable(model, f[1:A] >= 0.0)
    @variable(model, y[1:A], Bin)

    # objectif = flux sortant de s1 (poids de la plus petite classe)
    s1 = sources[1]
    @objective(model, Max, sum(f[a] for a in out_arcs[s1]))

    # ordre non décroissant des classes
    for i in 1:k-1
        si, sip1 = sources[i], sources[i+1]
        @constraint(model, sum(f[a] for a in out_arcs[si]) <= sum(f[a] for a in out_arcs[sip1]))
    end

    # conservation + consommation
    for v in 1:n
        @constraint(model, sum(f[a] for a in in_arcs[v]) - sum(f[a] for a in out_arcs[v]) == W[v])
    end

    # f <= M y, M = w(G)
    for a in 1:A
        @constraint(model, f[a] <= wG * y[a])
    end

    # au plus 1 arc depuis chaque source
    for si in sources
        @constraint(model, sum(y[a] for a in out_arcs[si]) <= 1)
    end

    # au plus 1 arc entrant par sommet
    for v in 1:n
        @constraint(model, sum(y[a] for a in in_arcs[v]) <= 1)
    end

    optimize!(model)

    obj  = objective_value(model)
    yval = value.(y)

    # reconstruction des classes à partir de y
    parent = fill(0, n)
    for v in 1:n
        for a in in_arcs[v]
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
    end

    return obj, classes
end

# ==========================
# 2) Formulation CUT-BASED (Section 2) + séparation (lazy constraints)
# ==========================
# Test de connectivité u~v après retrait d’un ensemble S
function connected_without_set(E::Array{Int,2}, u::Int, v::Int, Sset::BitVector)
    if Sset[u] || Sset[v]
        return false
    end
    n = size(E,1)
    seen = falses(n)
    q = Int[]
    push!(q, u); seen[u] = true
    while !isempty(q)
        x = popfirst!(q)
        if x == v
            return true
        end
        for y in 1:n
            if E[x,y] == 1 && !seen[y] && !Sset[y]
                seen[y] = true
                push!(q, y)
            end
        end
    end
    return false
end

# Séparateur minimal (brute force) — OK pour petites instances.
function min_vertex_separator_bruteforce(E::Array{Int,2}, u::Int, v::Int, cap::Vector{Float64})
    n = size(E,1)
    verts = [w for w in 1:n if w != u && w != v]
    bestS   = Int[]
    bestVal = Inf
    m = length(verts)
    for mask in 0:(UInt(1)<<m)-UInt(1)
        Sset = falses(n)
        val  = 0.0
        for j in 1:m
            if (mask >> (j-1)) & 0x1 == 0x1
                z = verts[j]
                Sset[z] = true
                val += cap[z]
                if val >= bestVal
                    break
                end
            end
        end
        if val >= bestVal
            continue
        end
        if !connected_without_set(E, u, v, Sset)
            bestVal = val
            empty!(bestS)
            for z in 1:n
                if Sset[z]; push!(bestS, z); end
            end
        end
    end
    return bestS, bestVal
end

function solve_bcpk_cut(E::Array{Int,2}, W::Vector{Int}, k::Int; tol::Float64=1e-6)
    n = length(W)

    # paires non adjacentes (u < v)
    nonedges = Tuple{Int,Int}[]
    for u in 1:n, v in u+1:n
        if E[u,v] == 0
            push!(nonedges, (u,v))
        end
    end

    model = Model(Gurobi.Optimizer)
    set_silent(model)

    # x[v,i] ∈ {0,1}
    @variable(model, x[1:n, 1:k], Bin)

    # objectif : poids de la classe 1 (et ordre forcé)
    @objective(model, Max, sum(W[v] * x[v,1] for v in 1:n))

    # ordre non décroissant des poids de classes
    for i in 1:k-1
        @constraint(model, sum(W[v]*x[v,i]   for v in 1:n) <=
                            sum(W[v]*x[v,i+1] for v in 1:n))
    end

    # chaque sommet dans au plus une classe
    for v in 1:n
        @constraint(model, sum(x[v,i] for i in 1:k) <= 1)
    end

    # Callback de séparation
    function conn_callback(cb_data)
        # valeurs courantes (fractionnelles possibles)
        xval = [JuMP.callback_value(cb_data, x[v,i]) for v in 1:n, i in 1:k]

        for i in 1:k
            cap = [xval[v,i] for v in 1:n]
            for (u,v) in nonedges
                s = xval[u,i] + xval[v,i]
                if s <= 1.0 + tol
                    continue
                end
                Smin, cutval = min_vertex_separator_bruteforce(E, u, v, cap)
                if s - cutval > 1.0 + tol
                    con = @build_constraint(x[u,i] + x[v,i] - sum(x[z,i] for z in Smin) <= 1)
                    MOI.submit(model, MOI.LazyConstraint(cb_data), con)
                end
            end
        end
        return
    end

    MOI.set(model, MOI.LazyConstraintCallback(), conn_callback)

    optimize!(model)

    obj  = objective_value(model)
    xopt = value.(x)

    classes = [Int[] for _ in 1:k]
    for i in 1:k, v in 1:n
        if xopt[v,i] > 0.5
            push!(classes[i], v)
        end
    end

    return obj, classes
end

# ==========================
# Démo avec l’instance fournie
# ==========================
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
# Symétriser (l’instance est non orientée)
E = max.(E, E')

k = 2

println("=== Méthode 1 : FLOW (Gurobi) ===")
objF, classesF = solve_bcpk_flow(E, W, k)
println("Poids minimal maximisé = ", objF)
for i in 1:k
    println("Classe $i : ", sort(classesF[i]))
end

println("\n=== Méthode 2 : CUT-BASED + séparation (Gurobi) ===")
objC, classesC = solve_bcpk_cut(E, W, k)
println("Poids minimal maximisé = ", objC)
for i in 1:k
    println("Classe $i : ", sort(classesC[i]))
end
