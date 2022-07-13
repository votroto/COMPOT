using JuMP
using MosekTools
using Combinatorics
using LinearAlgebra
using Gurobi

function _default_optimizer()
        attributes = [
                "QUIET" => true
        ]
        optimizer_with_attributes(Mosek.Optimizer, attributes...)
end

function _nonconvex_gurobi()
        attributes = [
		"NonConvex" => 2,
        ]
        optimizer_with_attributes(Gurobi.Optimizer, attributes...)
end

function edges(c)
	lc = length(c)
	rr = [[(c[i], c[i + 1]) for i in 1:(lc-1)]; (c[lc],c[1])]

end

function cycle(xs, n)
	@show cs = collect(combinations(xs, n))
	@show collect(map(edges, cs))
end

function generate_matrix(n; optimizer=_default_optimizer())
	combs = combinations

        N = 1:n
        W = rand(0:1, n, n)
        m = Model(optimizer)

        @variable m 0 <= P[N, N] <= 1
        @objective m Min dot(P, W)
        @constraints(m, begin
		[c in cycle(N, 3)], 1 <= sum(P[i,j] for (i,j) in c) <= 2
		[c in cycle(N, 3)], 1 <= sum(P[i,j] for (j,i) in c) <= 2
        	[(i, j) in combs(N, 2)], P[i, j] + P[j, i] == 1
        	[i in N], P[i, i] == 0
	end)

        optimize!(m)

        collect(value.(P)), termination_status(m)
end

function solve_feasibility(P; optimizer=_default_optimizer())
	combs = combinations

        candidates = axes(P, 1)
        perms = collect(permutations(candidates))

        m = Model(optimizer)
        @variable m w[perms] >= 0
        @constraints(m, begin
		[π in perms, (i, j) in combs(π, 2)], w[π] == P[i, j]
	end)

        optimize!(m)

        value.(w), termination_status(m)
end

function solve_qp(n; optimizer=_nonconvex_gurobi())
	combs = combinations

        N = 1:n
        perms = collect(permutations(N))

        m = Model(optimizer)

        @variable m 0 <= P[N, N]
	@variable m X[N, N] <= 1
        @objective m Max dot(P,X)
        @constraints(m, begin
        	[i in N], P[i, i] == 0
        	[(i, j) in combs(N, 2)], P[i, j] + P[j, i] == 1
		[c in cycle(N, 3)], 1 <= sum(P[i,j] for (i,j) in c) <= 2
		[π in perms], sum(X[i,j] for (i, j) in combs(π, 2)) <= 0
	end)

	optimize!(m)

	p_val = collect(value.(P))
	x_val = collect(value.(X))
	objective_value(m), termination_status(m), p_val, x_val
end


function one_sample()
	@show mat, status_mat = generate_matrix(3)
	@show w, status_feas = solve_feasibility(mat)
end