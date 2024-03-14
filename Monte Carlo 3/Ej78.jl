using Plots, Random
using ShiftedArrays
using Statistics, StatsBase

function gauss(y , delta)
    x = y + delta * (2 * rand() - 1) 

    # Accept x with probability  min(1, q(x|y))
    if abs(x) > abs(y)
        if(rand() > exp(0.5*(y-x)*(y + x)))
            return y # Rejection
        end
    end

    return x # Accept
end

function corrFunc(gi)
    kArr = [i for i in 1:1000]
    ρ = []
    meanVal = mean(gi)
    sigma = (sum(gi.^2) / N - (sum(gi) / N)^2 )
    for k in kArr
        # Prepare the shifted series
        giShift = ShiftedArrays.lead(gi, k)
        giShift = filter(!ismissing, giShift)
        M = length(giShift)
        giComp = gi[1:M]

        value = mean(giComp .* giShift) - meanVal * mean(giShift)
        push!(ρ, value/sigma)
    end

    return kArr, ρ
end

function gaussiana(x, μ, σ)
    return exp(-(x - μ)^2 / (2 * σ^2)) / sqrt(2 * π * σ^2)
end

Random.seed!(3688889)

x0 = (rand() - 1) # Uniform number between 1, -1
N = 10^6
delta = 0.1

xrand = []

for i in 1:N
    x0 = gauss(x0, delta)
    push!(xrand, x0)
end

x = range(-5, stop=5, length=100)
y = gaussiana.(x, 0, 1)

histogram(xrand, normalize = true, bins = 100)
plot!(x, y, label="Gaussiana", xlabel="x", ylabel="f(x)", title="Función Gaussiana")

g1 = xrand .^2
g2 = cos.(xrand)

k1, corr1 = corrFunc(g1)
k2, corr2 = corrFunc(g2)

plot(k1, corr1, label = "\$ G_1(x) = x^2 \$", xlabel = "\$ k\$", ylabel = "\$ ρ_G(k)\$")
plot!(k2, corr2, label = "\$ G_2(x) = cos(x) \$")

root = joinpath(@__DIR__, "rho.png")
savefig(root)

# N = 10^7
# delta = 0.1
tau1= sum(corr1) # 306.400
tau2 = sum(corr2) # 291.336

tau1Approx = corr1[1] / (1 - corr1[1]) # 311.317
tau2Approx = corr2[1] / (1 - corr2[1]) # 285.625

# ---------- Computing the integrals ----------
I1 = mean(g1)
I2 = mean(g2)

sigma1 = (sqrt(sum(g1.^2) / N - (sum(g1) / N)^2 ) * sqrt(2*tau1 + 1)) / sqrt(N)
sigma2 = (sqrt(sum(g2.^2) / N - (sum(g2) / N)^2 ) * sqrt(2*tau2 + 1)) / sqrt(N)

sigma1noC = sqrt(sum(g1.^2) / N - (sum(g1) / N)^2 )  / sqrt(N)
sigma2noC = sqrt(sum(g2.^2) / N - (sum(g2) / N)^2 )  / sqrt(N)
print(I1," ; ", sigma1)
print(I2," ; ", sigma2)