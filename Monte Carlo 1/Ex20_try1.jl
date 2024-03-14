using Distributions, Plots, Statistics, Random, HCubature

#Function to integrate
fun(x::Real) = sqrt(x) * cos(x)

#Function for the Uniform Sampling
function imp_sampling(M)
    G_list = []
    for i in 1:M
        x = -log(rand(Uniform(0, 1)))
        push!(G_list, fun(x) )
    end

    int = sum(G_list) / M
    sigma = sqrt( sum(G_list.^2) / M - (sum(G_list) / M)^2 ) / sqrt(M)
    
    return int, sigma
end

M = 10^6
imp_sampling(M)