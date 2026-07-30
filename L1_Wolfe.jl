start=time()

using Pkg
Pkg.instantiate()

using Distributed
using SlurmClusterManager   # separate package

# Start workers (blocks until started)
addprocs(SlurmManager())

# if nworkers() == 1
#     addprocs(44, exeflags="--project=$(Base.active_project())")
# end

@everywhere using JuMP, JLD, LinearAlgebra, Polyhedra, Random, Distributions, CDDLib, Base.Threads

@everywhere BLAS.set_num_threads(1)

@sync @everywhere begin
    pid = getpid()
    # current CPU (Linux): returns the CPU this thread is currently running on
    cpu = try
        ccall(:sched_getcpu, Int32, ())
    catch
        -1
    end
    println("worker myid=$(myid()) pid=$pid cpu=$cpu")
end

println("Done setup in time ", time() - start)

@everywhere include("exact_frechet_computation.jl")

args = parse.(Int, ARGS)

@everywhere k_range =collect($args[1]:$args[2])
@everywhere n_range =collect($args[3]:$args[4])
@everywhere reps =$args[5]

# # Pre-compute all alphas before threading
@everywhere alphas_dict = Dict(k => precompute_l1_alphas(k) for k in k_range)

@everywhere total_iterations = length(k_range) * length(n_range) * reps

println("Defined functions in time ", time()-start)

pmap(x -> Wolfe_experiment_runner(x, k_range, n_range, reps, alphas_dict), 1:total_iterations)

println("All done in time ", time()-start)

dimensions = -2*ones(Int,(reps, length(n_range), length(k_range)))

# for result in results
#     rep, n_id, k_id, dim = result
#     dimensions[rep, n_id, k_id] = dim
# end

# println("Dimensions array:", dimensions)

# JLD.save("l1.jld", "dimensions", dimensions)
println("Computation complete. Results saved to l1.jld")