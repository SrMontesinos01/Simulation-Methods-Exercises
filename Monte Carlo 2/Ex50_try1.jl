using Plots, SpecialFunctions, Distributions, HCubature

# Original Probability Desity Function
f_n(x::Float64) = exp(+0.5 * x^2 - x^4)
f(x::Float64) = exp(+0.5 * x^2 - x^4) / C
hg(x::Float64) = C2* exp(+0.5 * x^2 - x^4) / sqrt(2π)
# Probability Function
h(x::Float64) = exp(x^2 - x^4 - 1/4)

h(0.1)
function rejectionRep()
    x = 0
  
    while true
        x = rand(Normal(0, 1))
        (rand() < h(x)) && break
    end
    
    return x
end

C2 = ℯ^(-1/4)
# Computing the normalization constant
C, err = hquadrature(f_n, -300, 300)
a = 1 / C
err_C = (1/C)^2 * err
round.(a, digits = 10)
round.(err_C , digits = 12)

sol  = hquadrature(hg, -300, 300)
acc, err = round.(sol, digits = 12)
C_acc =  C2/ (acc * sqrt(2π))
err_acp = (1/acc)^2 * err / sqrt(2π)
round.(C_acc, digits = 10)
round.(err_acp , digits = 12)

# Obtaining the random numbers using the rejection method
M = 10^6
x = collect(-2:0.1:2)

a = [rejectionRep() for i in 1:M]

histogram(a, bins = 100, normalize=:pdf, color = "lightgreen",
            title="M = $M", label = "Histogram", xlabel = "x", ylabel ="Probability")
plot!(x, f.(x),  linewidth=2, linestyle=:dash, color = "red", label = "\$ f_{̂x}(x) \$")

root = joinpath(@__DIR__, "images")
name = "\\Ex50_$M.png"
savefig(root * name)