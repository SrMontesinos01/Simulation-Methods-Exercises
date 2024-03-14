using Plots, SpecialFunctions, Distributions

# Original Probability Desity Function
f(x::Float64, α::Float64) = (x^(α - 1) * exp(-x)) / gamma(α)
# New Probability Desity Function
g(x::Float64, α::Float64) = x <= 1 ? x^(α - 1) / gamma(α) :  exp(-x) / gamma(α)
# Probability Function
h(x::Float64, α::Float64) = x <= 1 ? exp(-x) : x^(α - 1)  

# Function that generates random numbers according to the pdf g(x)
function rand_g(α::Float64, limit::Float64, C::Float64)
    u = rand()

    if u <= limit
        return (u * gamma(α + 1) *C)^(1/α)
    else
        return -log(-u * C * gamma(α) + 1/α + 1/ℯ)
    end
end

function rejectionRep_47(alpha::Float64, limit::Float64, C::Float64)
    x = 0
    while true
        x = rand_g(α, limit, C)
        u = rand()
        u < h(x, alpha) && break
    end
    
    return x
end

# Obtaining the random numbers using the rejection method
M = 10^6
α = 0.4
x = collect(0:0.0085:4)

C = (α + ℯ) / (gamma(α + 1) * ℯ)
limit = 1 /(gamma(α + 1) * C) # Cumulative at x = 1

A = [rejectionRep_47(α, limit, C) for i in 1:M]
histogram(A, bins = 300, normalize=:pdf, color = "lightgreen", xlim = (0,3.0),
            title="α= $α   ;  M = $M", label = "Histogram", xlabel = "x", ylabel ="Probability")
plot!(x, f.(x,α),  linewidth=2, linestyle=:dash, color = "red", label = "\$ f_{̂x}(x) \$")

root = joinpath(@__DIR__, "images")
name = "\\Ex47_$α.png"
savefig(root * name)