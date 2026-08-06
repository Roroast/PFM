start=time()

using Pkg
Pkg.instantiate()

using Distributed
using SlurmClusterManager   # separate package

# Start workers (blocks until started)
addprocs(SlurmManager())

include("exact_frechet_computation.jl")
args = parse.(Int, ARGS)
k_range = collect(args[1]:args[2])
n_range = collect(args[3]:args[4])
reps = args[5]
alphas_dict = Dict(k => precompute_l1_alphas(k) for k in k_range)

function setup_workers(ids)
    @everywhere ids begin
        using JuMP, JLD, LinearAlgebra, Polyhedra, Random, Distributions, CDDLib, Base.Threads
        BLAS.set_num_threads(1)
        include("exact_frechet_computation.jl")
    end
    @everywhere ids k_range = $k_range
    @everywhere ids n_range = $n_range
    @everywhere ids reps = $reps
    @everywhere ids alphas_dict = Dict(k => precompute_l1_alphas(k) for k in k_range)
end

println("Initial setup...")
setup_workers(workers())

total_iterations = length(k_range) * length(n_range) * reps

println("Defined functions and setup workers in time ", time() - start)

task_queue = Channel{Int}(total_iterations + 1000)
for i in 1:total_iterations; put!(task_queue, i); end

completed = 0
c_lock = ReentrantLock()
num_initial_workers = nworkers()

@sync begin
    for w_initial in workers()
        @async begin
            w = w_initial
            while true
                t = nothing
                lock(c_lock) do
                    if completed < total_iterations && isready(task_queue)
                        t = take!(task_queue)
                    end
                end
                t === nothing && break
                
                try
                    remotecall_fetch(Wolfe_experiment_runner, w, t, k_range, n_range, reps, alphas_dict)
                    lock(c_lock) do
                        completed += 1
                        if completed % 10 == 0
                            println("Progress: $completed / $total_iterations tasks completed.")
                        end
                        if completed % 100 == 0
                            GC.gc() # Periodically clean master
                        end
                    end
                catch e
                    if e isa ProcessExitedException
                        println("Worker $w died during task $t. Replacing...")
                        put!(task_queue, t)
                        new_ws = addprocs(SlurmManager())
                        if !isempty(new_ws)
                            setup_workers(new_ws)
                            w = new_ws[1]
                        else
                            println("Could not replace worker $w. This slot will remain empty.")
                            break
                        end
                    else
                        println("Error on worker $w for task $t: $e")
                        put!(task_queue, t)
                        rethrow(e)
                    end
                end
            end
        end
    end
end

println("All done in time ", time()-start)

dimensions = -2*ones(Int,(reps, length(n_range), length(k_range)))

# for result in results
#     rep, n_id, k_id, dim = result
#     dimensions[rep, n_id, k_id] = dim
# end

# println("Dimensions array:", dimensions)

# JLD.save("l1.jld", "dimensions", dimensions)
println("Computation complete. Results saved to l1.jld")