using Plots
using LaTeXTabulars
using LaTeXStrings


function dimensions_from_files(prob::String, n_range::Vector{Int}, k_range::Vector{Int}, rep_count::Int)
    
    max_n = n_range[end]
    max_k = k_range[end]
    folder_name = prob*"_$(max_n)_$(max_k)_$(rep_count)"
    file_names = readdir(folder_name)

    err_files = file_names[findall(endswith.(file_names, ".err"))]
    out_files = file_names[findall(endswith.(file_names, ".out"))]

    for file in err_files
        if read(folder_name*"/"*file, String) != ""
            print(file)
        end
    end

    dimensions = -2*ones(Int,(rep_count, length(n_range), length(k_range)))
    timings = -2*ones(Float64,(rep_count, length(n_range), length(k_range)))
    convergeregex = r"k = (\d+), n = (\d+), rep = (\d+), dim = (\d+), time = (\d+\.\d+) seconds"

    for file in out_files
        str = split(read(folder_name*"/"*file, String), "\n")
        str = str[findall(startswith.(str, "k = "))]
        for result in str
            data = match(convergeregex, result).captures[:]
            k, n, r, d = parse.(Int,data[1:4])
            t = parse(Float64,data[5])
            dimensions[r,findfirst(n_range .== n),findfirst(k_range .== k)] = d
            timings[r,findfirst(n_range .== n),findfirst(k_range .== k)] = t
        end
    end

    error_experiments = findall([-2*ones(Int64,length(n_range),length(k_range)) == dimensions[i,:,:] for i in 1:rep_count])

    average_timings = ((sum(timings[:, :, :], dims=1)/rep_count)[1,:,:])/60

    latex_tabular("LatexTables/"*prob*"_average_times.tex",
        Tabular("l|"*"l"^(length(k_range))),
        [Rule(:top),
        [0, k_range...],
        Rule(),
        hcat("n=".*string.(n_range),round.(average_timings,sigdigits=2)),
        Rule(:bottom)]
    )

    return error_experiments, dimensions, timings
end

function plotting(dimensions, problem, plot_type::String, n_range::Vector{Int}, k_range::Vector{Int}, rep_count::Int, 
    size = (800,600),
    error_experiments::Vector{Int}=[0],
    leg = :topleft)

    dimensions = dimensions[setdiff(1:rep_count, error_experiments),:,:]
    rep_count = rep_count - length(error_experiments)

    if plot_type == "ex_dim"
        data = transpose(transpose(sum(dimensions[:, :, :], dims=1)[1,:,:])./((k_range.-1)*rep_count))
    elseif plot_type == "uniqueness_prob"
        data = (sum(dimensions[:, :, :] .== 0, dims=1)/rep_count)[1,:,:]
    end

    P = plot(xlabel=L"Sample Size, $n$", size = size)

    for (k_id, k) in enumerate(k_range)
        plot!(n_range, vec(data[:, k_id]), 
        xticks = n_range,
        label = string("\$k=",k,"\$"),
        bottommargin = 5Plots.mm,
        leftmargin = 10Plots.mm,
        leg = leg) 
    end

    if plot_type == "ex_dim"
        if problem == "L_1"
            title!(L"Average Dimension of $\ell_1$ Fréchet Mean Set")
        elseif problem == "Trop"
            title!("Average Dimension of Tropical Fréchet Mean Set")
        elseif problem == "L_{inf}"
            title!(L"Average Dimension of $\ell_{\infty}$ Fréchet Mean Set")
        end
        ylabel!(L"Average Fréchet Mean Set Dimension $/(k-1)$")
        savefig(P,"Figures/"*problem*"_expected_dimension.pdf")
    elseif plot_type == "uniqueness_prob"
        if problem == "L_1"
            title!(L"Empirical Probability of a Unique $\ell_1$ Fréchet Mean")
        elseif problem == "Trop"
            title!("Empirical Probability of a Unique Tropical Fréchet Mean")
        elseif problem == "L_{inf}"
            title!(L"Empirical Probability of a Unique $\ell_{\infty}$ Fréchet Mean")
        end
        ylabel!(L"Empirical Uniqueness Probability, $p_n$")
        savefig(P,"Figures/"*problem*"_uniqueness_probability.pdf")
    end

    latex_tabular("LatexTables/"*problem*"_"*plot_type*".tex",
        Tabular("l|"*"l"^(length(k_range))),
        [Rule(:top),
        [0, k_range...],
        Rule(),
        hcat("n=".*string.(n_range),round.(data,sigdigits=2)),
        Rule(:bottom)]
    )

    return data

end