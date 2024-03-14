using Random, DelimitedFiles, Plots
using Statistics, StatsBase
using  ShiftedArrays

# Seed for generating Random numbers
Random.seed!(34126)

# Defining the functions of the SDE: ̇x(t) = -ax + √Df(x)ξ_{OU}(t)
q(x) = - 2*x
g(x, D::Float64) = sqrt(D) .* f(x)
f(x) = 1

# Function for integrating the SDE
function heunOrnsteinUhlenbeck(x, t, h , τ, nwrt, nstp; name = "tray.txt")
    """
    Function that performs the Heun method for Ornstein-
    Uhlenbeck noise. 

    Given a vector of Initial conditions x, write
    the trajectory for each initial condition in
    "tray.txt" (or corresponding name)

    """
    # Computing auxiliar values
    n = length(x)
    hs = sqrt(h)
    p = exp(-h/τ)
    beta = - τ * (1 - p) / hs
    gamma = sqrt(0.5*τ*(1-p^2) - beta^2)
    aOld, bOld, ghOld = zeros(n), sqrt(0.5*τ)*(p -1) .* randn(n), zeros(n)

    root = joinpath(@__DIR__, name)
    file = open(root, "w")
    for i in 1:nwrt
        #Write results
        write(file, "$t")
        for num in x
            write(file, "\t$num")
        end
        write(file, "\n")

        for j in 1:nstp
            u = randn(n)
            a, b = hs.*u, beta.*u + gamma*randn(n)
            gh = p.*(ghOld - aOld) + a - bOld + b
            aOld, bOld, ghOld = a, b, gh

            k = h.*q.(x) + gh .* g(x, D)
            x = x .+ 0.5*(k .+ h*q.(x + k) .+ gh .*g(x + k, D))
            t += h
        end
    end
    close(file)
end

# Function for b)
function tauPlot(tauData, tauList)
    plots = []
    for (i, data) in enumerate(tauData)

        xcols = size(data)[2] - 1 # Number of columns with x(t) data
        tau = tauList[i]
        p = plot(xlabel = "\$t\$", ylabel = "\$ x(t) \$", label = false) # plot for a constant τ
        title!("\$τ =$tau\$", fontsize=0.5)
        print("holalaa:", xcols)
        for col in 2:xcols
            print(col)
            
            p = plot!(data[:,1], data[:,col], label = "",dpi = 1000)
        end
        push!(plots, p)
    end

    pf = plot(plots..., layout = (2,2))

    return pf
end

# Function for c),d)
function calcMoment(nData, nList; exp = 1)
    nAve = []
    for (i, data) in enumerate(nData)
        n = nList[i]
        rows = size(data)[1]
        xAve= []
        print(size(data))
        for j in 1:rows
            μ = mean(data[j, 2:end].^exp)
            push!(xAve, μ)
        end 

        push!(nAve, xAve)
    end

    return nAve
end

# Function for e), f)
function calcCorr(data, timeStep::Int64)
    """
    data should be a matrix-like object with each time series as a column
    """
    # 20 = 1/(h * nstp)
    sArr = [i for i in 0.05:0.05:5] * 20 # Lags for the correlation
    corr = []
    for s in sArr
        # Prepare the shifted series
        # xt_shift = ShiftedArrays.lead(xt, Int64(trunc(5)), default = 0)
        xt = data[timeStep, 2:end]
        timeStepShifted = Int64(trunc(timeStep + s))
        xtShift = xtData[timeStepShifted, 2:end]

        push!(corr ,mean(xt .* xtShift))
    end

    return sArr/20, corr
end

# Function for e), f)
function corrTimePlot(data, timeArr; tau = 0.5)
    corrTimes = []
    p = plot(xlabel = "\$s\$", ylabel = "\$C_n(t,s)\$", legendfont = 10)
    title!("\$τ =$tau\$", fontsize=0.5)
    for time in timeArr
        timeStep = time / (h * nstp)
        sArr, corrAve = calcCorr(data, Int64(trunc(timeStep)))
        corrAve = corrAve ./corrAve[1]
        p = plot!(sArr, corrAve, label = "\$ C_n(t = $time,s) \$", dpi = 1000)

        # Obtains the correlation time for each "time"
        indexCorrTime = argmin(abs.(corrAve .- exp(-1)))
        print(indexCorrTime, "\n")
        corrTime = sArr[indexCorrTime]
        print(corrTime, "\n")
        push!(corrTimes, corrTime)
    end

    return p, corrTimes
end

# ---------------------------------------------------- #
# --------------------     b)    --------------------- #
# ---------------------------------------------------- #
# Defining parameters
h = 0.01
nwrt, nstp = 50, 10
totalTime = h*nwrt*nstp
D = 0.01

t0, x0 = 0.0, [1.0 for i in 1:10] # t0, x0
tauList = [2, 0.5, 0.05, 0.001]

# Simulation for several values of τ
for (i, tau) in enumerate(tauList)
    fname = "tray_τ$i.txt"
    heunOrnsteinUhlenbeck(x0, t0, h , tau, nwrt, nstp, name = fname)
end

# Reading the data
xtDataTau = []
for (i, tau) in enumerate(tauList)
    root = joinpath(@__DIR__, "tray_τ$i.txt")
    xtData = readdlm(root, '\t', Float64) # reading all the trayectories for a certain tau
    push!(xtDataTau, xtData)
end

# Performing the plot
p = tauPlot(xtDataTau, tauList)
p = plot!(size = (600, 600))
root = joinpath(@__DIR__, "tray_b.png")
savefig(root)

# ---------------------------------------------------- #
# --------------------     c)    --------------------- #
# ---------------------------------------------------- #
nTray = [10, 100, 1000]
tau = 0.5
# Simulation for several amount of trajectories
for n in nTray
    t0, x0 = 0.0, [1.0 for i in 1:n] # t0, x0
    fname = "tray_n$n.txt"
    heunOrnsteinUhlenbeck(x0, t0, h , tau, nwrt, nstp, name = fname)
end

# Reading the data
nTrayData = []
for (i, n) in enumerate(nTray)
    root = joinpath(@__DIR__, "tray_n$n.txt")
    xtData = readdlm(root, '\t', Float64) # reading all the trayectories for a certain tau
    push!(nTrayData, xtData)
end

nAve = calcMoment(nTrayData, nTray)
t = nTrayData[1][:,1]

p1 = plot(t ,nAve[1], label = "\$ n = 10\$", ylabel = "\$ ⟨x(t)⟩\$" , xlabel = "\$ t \$")
p1 = plot!(t ,nAve[2], label = "\$ n = 100\$")
p1 = plot!(t ,nAve[3], label = "\$ n = 1000\$")

p2 = plot(t ,nAve[1], label = "\$ n = 10\$", ylabel = "\$ ⟨x(t)⟩\$" , xlabel = "\$ t \$")
p2 = plot!(t ,nAve[2], label = "\$ n = 100\$")
p2 = plot!(t ,nAve[3], label = "\$ n = 1000\$")
xlims!(1.0, 5)
ylims!(-0.01, 0.05)

pf = plot(p1, p2, layout = (1,2))
pf = plot!(size = (600, 300))
root = joinpath(@__DIR__, "tray_c.png")
savefig(root)
# ---------------------------------------------------- #
# --------------------     d)    --------------------- #
# ---------------------------------------------------- #
xDet(t::Float64) = exp(-2*t) # Deterministic Solution

nTray = [10, 100, 1000]
tau = 0.5
# Simulation for several amount of trajectories
for n in nTray
    t0, x0 = 0.0, [1.0 for i in 1:n] # t0, x0
    fname = "tray_n$n.txt"
    heunOrnsteinUhlenbeck(x0, t0, h , tau, nwrt, nstp, name = fname)
end

# Reading the data
nTrayData = []
for (i, n) in enumerate(nTray)
    root = joinpath(@__DIR__, "tray_n$n.txt")
    xtData = readdlm(root, '\t', Float64) # reading all the trayectories for a certain tau
    push!(nTrayData, xtData)
end

nAve = calcMoment(nTrayData, nTray, exp = 2)
t = nTrayData[1][:,1]
x2 = xDet.(t) .^2
# Plotting results
p1 = plot(t ,nAve[1], label = "\$ n = 10\$", ylabel = "\$ ⟨x^2(t)⟩\$" , xlabel = "\$ t \$")
p1 = plot!(t ,nAve[2], label = "\$ n = 100\$")
p1 = plot!(t ,nAve[3], label = "\$ n = 1000\$")
p1 = scatter!(t ,x2 , label = "\$ x^2_{det}(t)\$", linestyle=:dash, 
                color =:orange, markersize = 2.0, markerstrokewidth = 0.5)

p2 = plot(t ,nAve[1], label = "\$ n = 10\$", ylabel = "\$ ⟨x^2(t)⟩\$" , xlabel = "\$ t \$")
p2 = plot!(t ,nAve[2], label = "\$ n = 100\$")
p2 = plot!(t ,nAve[3], label = "\$ n = 1000\$")
p2 = scatter!(t ,x2 , label = "\$ x^2_{det}(t)\$", linestyle=:dash, 
                color =:orange, markersize = 2.0, markerstrokewidth = 0.5)

xlims!(1, 5)
ylims!(-0.0005, 0.005)

pf = plot(p1, p2, layout = (1,2))
pf = plot!(size = (600, 300))
root = joinpath(@__DIR__, "tray_d.png")
savefig(root)
# ---------------------------------------------------- #
# --------------------     e)    --------------------- #
# ---------------------------------------------------- #
# Defining parameters
h = 0.01
nwrt, nstp = 250, 5
totalTime = h*nwrt*nstp

tau = 0.5
t0, x0 = 0.0, [1.0 for i in 1:1000] # t0, x0
fname = "tray_e.txt"
heunOrnsteinUhlenbeck(x0, t0, h , tau, nwrt, nstp, name = fname)

# Reading the data
root = joinpath(@__DIR__, "tray_e.txt")
xtData = readdlm(root, '\t', Float64) 

timeArr = [0.5, 1.0, 2.0, 5.0] # Stationary regime for time > 2 approx.
pf, corrTimes = corrTimePlot(xtData, timeArr)
pf
corrTimes

root = joinpath(@__DIR__, "tray_e.png")
savefig(root)
# ---------------------------------------------------- #
# --------------------     f)    --------------------- #
# ---------------------------------------------------- #
# Defining parameters
h = 0.01
nwrt, nstp = 250, 5
totalTime = h*nwrt*nstp
tauList = [0.05, 2.0]

t0, x0 = 0.0, [1.0 for i in 1:1000] # t0, x0
# Simulation for several values of τ
for (i, tau) in enumerate(tauList)
    fname = "tray_e_τ$i.txt"
    heunOrnsteinUhlenbeck(x0, t0, h , tau, nwrt, nstp, name = fname)
end

# Reading the data
tauCorrTimes = []
tauPlots = []
timeArr = [0.5, 1.0, 2.0, 5.0]
for (i, tau) in enumerate(tauList)
    root = joinpath(@__DIR__, "tray_e_τ$i.txt")
    xtData = readdlm(root, '\t', Float64) # reading all the trayectories for a certain tau
    pf, corrTimes = corrTimePlot(xtData, timeArr, tau = tau)
    push!(tauCorrTimes, corrTimes)
    push!(tauPlots, pf)
end

tauCorrTimes
pf = plot(tauPlots..., layout = (1,2),  size=(750, 400))
pf = plot!(size = (700, 400))

root = joinpath(@__DIR__, "tray_f.png")
savefig(root)
# ---------------------------------------------------- #
# --------------------     g)    --------------------- #
# ---------------------------------------------------- #
f(x) = x

# Defining parameters
h = 0.01
nwrt, nstp = 500, 1
totalTime = h*nwrt*nstp
D = 0.01
tau = 0.05

t0, x0 = 0.0, [1.0 for i in 1:10] # t0, x0
tauList = [0.5]

# Simulation for several values of τ
for (i, tau) in enumerate(tauList)
    fname = "tray_g_τ$i.txt"
    heunOrnsteinUhlenbeck(x0, t0, h , tau, nwrt, nstp, name = fname)
end

# Reading the data
xtDataTau = []
for (i, tau) in enumerate(tauList)
    root = joinpath(@__DIR__, "tray_g_τ$i.txt")
    xtData = readdlm(root, '\t', Float64) # reading all the trayectories for a certain tau
    push!(xtDataTau, xtData)
end

xtDataTau[1]

# Performing the plot
p1 = tauPlot(xtDataTau, tauList)
p2 = tauPlot(xtDataTau, tauList)
xlims!(2.5, 5)
ylims!(-0.0005, 0.005)

pf = plot(p1, p2, layout = (1,2))
pf = plot!(size = (600, 300))

title!("")
root = joinpath(@__DIR__, "tray_g.png")
savefig(root)
