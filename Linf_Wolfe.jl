start=time()

using Pkg
Pkg.instantiate()

using Distributed
using SlurmClusterManager   # separate package

# Start workers (blocks until started)
addprocs(SlurmManager())

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
@everywhere r_range =collect($args[5]:$args[6])

# # Pre-compute all alphas before threading
@everywhere alphas_dict = Dict(k => precompute_linf_alphas(k) for k in k_range)

@everywhere total_iterations = length(k_range) * length(n_range) * length(r_range)

@everywhere before = Base.gc_num()

println("Defined functions in time ", time()-start)

pmap(x -> Wolfe_experiment_runner(x, k_range, n_range, r_range, alphas_dict), 1:total_iterations)

println("All done in time ", time()-start)