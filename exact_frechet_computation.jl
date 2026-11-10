using JuMP, LinearAlgebra, CDDLib, Polyhedra

function lifted_Q(data::Matrix{Rational{BigInt}},alphas::Matrix{Int64})

    k=size(data,1);
    n=size(data,2);
    m = size(alphas,1);

    A = vcat([hcat(vcat([hcat([i==l for i in 1:n]...) for j in 1:m]...),alphas) for l in 1:n]...);
    b = -vec(alphas*data);

    h = LiftedHRepresentation{Rational{BigInt}}(hcat(b,A))
    p = Polyhedra.polyhedron(h, CDDLib.Library(:exact))

    return p

end

function distances(data::Matrix{Rational{BigInt}}, theta::Vector{Rational{BigInt}}, alphas::Matrix{Int64})
    return maximum(alphas * (data .- theta),dims = 1)
end

function distances_with_indices(data::Matrix{Rational{BigInt}}, theta::Vector{Rational{BigInt}}, alphas::Matrix{Int64})
    dist = maximum(alphas * (data .- theta),dims = 1)
    return dist,[findall(dist[i] .== (alphas * (data .- theta))[:,i]) for i in 1:size(data)[2]]
end

function precompute_linf_alphas(k::Int64)
    return vcat(Matrix(+I,k,k),Matrix(-I,k,k));
end

function precompute_l1_alphas(k::Int64)
    return vcat(collect.(vcat(collect(Iterators.product(ntuple(i->[-1,1], k)...))...))'...);
end

function precompute_tropical_alphas(k::Int64)
    id_matrix = Matrix(I, k, k)
    pairwise = vcat([[(i == a) - (i == b) for i in 1:k]' for a in 1:k for b in 1:k if a != b]...)
    return vcat(id_matrix, -id_matrix, pairwise)
end

function Wolfe_FM_dimension(data::Matrix{Rational{BigInt}}, alphas::Matrix{Int64})
    
    p = lifted_Q(data, alphas);
    
    m = size(alphas,1);
    n = size(data,2);
    k = size(data,1);

    start_distances = distances(data, zeros(Rational{BigInt}, k), alphas)
          
    optimum_distances, theta = lifted_wolfe(p, n, k, start_distances);
    
    b = Rational.(vec(-maximum(alphas*data - hcat([vcat([optimum_distances[j] for i in 1:m]...) for j in 1:n]...),dims=2)));
    h = LiftedHRepresentation{Rational{BigInt}}(hcat(b,alphas))

    FMset = Polyhedra.polyhedron(h, CDDLib.Library(:exact))

    solution_checker(data,alphas,FMset,theta,optimum_distances)

    dimension = Polyhedra.dim(FMset)

    return dimension

end

function Wolfe_experiment_runner(flat_idx::Int64, k_range::Vector{Int64}, n_range::Vector{Int64}, r_range::Vector{Int64}, alphas_dict::Dict{Int64, Matrix{Int64}})
    
    # Convert flat index back to (k_id, n_id)
    r_id = div(flat_idx - 1, length(k_range)*length(n_range)) + 1
    n_id = div(rem(flat_idx - 1, length(k_range)*length(n_range)), length(k_range)) + 1
    k_id = rem(flat_idx - 1, length(k_range)) + 1
    
    rep = r_range[r_id]
    k = k_range[end - k_id + 1]
    n = n_range[end - n_id + 1]

    alphas = alphas_dict[k]
    
    #TODO: try other distributions, uniform, norm radial etc....
    data = Rational{BigInt}.(randn(k, n))

    start_time = time()

    dimension = Wolfe_FM_dimension(data, alphas)

    end_time = time()

    println("k = $k, n = $n, rep = $rep, dim = $dimension, time = $(round(end_time - start_time; digits=2)) seconds")

end

function solution_checker(data::Matrix{Rational{BigInt}}, alphas::Matrix{Int64}, FMset::CDDLib.Polyhedron{Rational{BigInt}}, theta::Vector{Rational{BigInt}}, optimum_distances::Vector{Rational{BigInt}})

    @assert in(theta,FMset)
    
    n = size(data,2)
    k = size(data,1)
    dim = n+k

    theta_distances, active_indices = distances_with_indices(data, theta, alphas)

    @assert theta_distances == optimum_distances'

    subgradients = [Polyhedra.polyhedron(convexhull([theta_distances[i]*alphas[j,:] for j in active_indices[i]]...), CDDLib.Library(:exact)) for i in 1:n]

    total_subgradient = sum(subgradients)

    @assert in(zeros(Rational{BigInt}, k),total_subgradient) "data was $data, FMset was $FMset"

end

function lifted_wolfe(p::CDDLib.Polyhedron{Rational{BigInt}},n::Int64,k::Int64,start_distances::Matrix{Rational{BigInt}})

    T = Rational{BigInt}
    total_dim = n+k

    model = GenericModel{T}(CDDLib.Optimizer{T})
    @variable(model, v[i=1:total_dim], start = vcat(start_distances', zeros(T, k))[i])
    @constraint(model, v in p)
    grad = zeros(T,total_dim)
    grad[1:n] = start_distances
    @objective(model, Min, sum(grad.*v))

    optimize!(model)
    p_vertex = value.(v)

    cap = total_dim+1
    C = hcat(p_vertex, zeros(T, total_dim, total_dim))
    m = 1
    lambda = vcat(1, zeros(T,cap))
    x_t = C[:,1:m] * lambda[1:m]

    ones_buf = ones(T,cap)
    G = zeros(T, 1, 1)

    grad[1:n] = x_t[1:n]
    set_objective_coefficient(model, v[1:n], grad[1:n])
    optimize!(model)
    p_vertex = value.(v)


    while dot(p_vertex[1:n],x_t[1:n]) < dot(x_t[1:n],x_t[1:n])

        m += 1
        C[:,m] = p_vertex
        lambda[m] = 0
        @assert all(lambda .>= 0)
        @assert (sum(lambda) == 1) "sum error 1"

        C1 = @view C[1:n, 1:m]
        A = zeros(T,m+1,m+1)
        G = zeros(T, m, m)
        y = zeros(T, m+1)
        mul!(G, transpose(C1), C1)
        b = zeros(T,m+1)
        b[1:m] .= 0
        b[m+1] = 1
        A[1:m,1:m] = G
	    A[1:m, m+1] = ones_buf[1:m]
        A[m+1, 1:m] = ones_buf[1:m]'

    	y = A\b
        y_coeffs = view(y, 1:m)
        @assert sum(y_coeffs) == 1 

        relevant_indices = y_coeffs .<= 0

        while any(relevant_indices)
	        idx = findall(relevant_indices)

            gamma,minimal_index = findmin(lambda[idx]./(lambda[idx] - y_coeffs[idx]))

            @. lambda[1:m] = (1-gamma) * lambda[1:m] + gamma * y_coeffs
            @assert all(lambda .>= 0)
            @assert sum(y_coeffs) == 1
            @assert (sum(lambda) == 1) "sum error 2"

    	    i = idx[minimal_index]
            @assert lambda[i] == 0 "lambda is $lambda, y_coeffs is $y_coeffs, gamma is $gamma, index is $minimal_index"

            C[:, i:(m-1)] = C[:,(i+1):m]
            C[:, m] = zeros(T,total_dim)
            lambda[i:(m-1)] = lambda[(i+1):m]
            lambda[m] = 0
            m -= 1

    	    C1 = @view C[1:n,1:m]
    	    A = zeros(T,m+1,m+1)
            G = zeros(T, m, m)
            y = zeros(T, m+1)
            mul!(G, transpose(C1), C1)
            A[1:m,1:m] = G
	        A[1:m, m+1] = ones_buf[1:m]
            A[m+1, 1:m] = ones_buf[1:m]'
            b = zeros(T,m+1)
            b[1:m] .= 0
            b[m+1] = 1

	        y = A\b
            y_coeffs = view(y, 1:m)
            @assert sum(y_coeffs) == 1 
            relevant_indices = y_coeffs .<= 0

        end

        lambda[1:m] = y_coeffs
        @assert all(lambda .>= 0)
        @assert (sum(lambda) == 1) "sum error 3"

        mul!(x_t, C[:,1:m], lambda[1:m])
	
	    grad[1:n] = x_t[1:n]
	    @assert all(grad[n+1:end] .== 0)

	    set_objective_coefficient(model, v[1:n], grad[1:n])
        optimize!(model)
        p_vertex = value.(v)

    end

    @assert (p_vertex[1:n]'*x_t[1:n]) >= (x_t[1:n]'*x_t[1:n])

    return x_t[1:n], x_t[n+1:end]

end