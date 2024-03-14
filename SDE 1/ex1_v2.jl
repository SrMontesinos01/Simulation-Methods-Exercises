module ex1_v2

using Distributions, Plots, Plots.PlotMeasures
using GLM, DataFrames, DelimitedFiles

# Function for implementing the milshtein algorithim
function milshtein(x, t, n_stp, h, h_sqrt, D, n)
    for j in 1:n_stp
        uh = h_sqrt * randn(n)
        uh2 = uh .^2
        x = x + h*q(x,t) + g(x,t, D)*(uh + 0.5 * gprime(x,t)*uh2)
        t = t + h
    end

    return x, t
end

# Function for b)
function moments_plot(xt_data)
    lim = size(xt_data)[1] # Dimension of each trayectory vector

    t = xt_data[:,1]
    x_ave = Vector{Float64}()
    x2_ave = Vector{Float64}()
    x4_ave = Vector{Float64}()
    σ2 = Vector{Float64}()
    for i in 1:lim 
        push!(x_ave, sum(xt_data[i,2:end]) / n)
        push!(x2_ave, sum(xt_data[i,2:end].^2) / n)
        push!(x4_ave, sum(xt_data[i,2:end].^4) / n)
        push!(σ2, x4_ave[i] - x2_ave[i]^2)
    end

    p1 =plot(t, x_ave, ylabel = "\$⟨x(t)⟩\$", xlabel = "t", label = "\$⟨x(t)⟩\$",
    color = "blue", legend=:topright, title = "\$⟨x(t)⟩\$")
    p2 = plot(t, x2_ave, ylabel = "\$⟨x^2(t)⟩\$", xlabel = "t",label = "\$⟨x^2(t)⟩\$",
    color = "green",legend=:bottomright, title = "\$⟨x^2(t)⟩\$")
    p3 = plot(t, x4_ave, ylabel = "\$⟨x^4(t)⟩\$", xlabel = "t", label = "\$⟨x^4(t)⟩\$",
    color = "red", legend=:bottomright, title = "\$⟨x^4(t)⟩\$")

    p4 =plot(t, x_ave, ylabel = "x(t) moments", xlabel = "t", label = "\$⟨x(t)⟩\$", 
    color = "blue", legend = false, title = "\$⟨x(t)⟩\$, \$⟨x^2(t)⟩\$, \$⟨x^4(t)⟩\$")
    plot!(t, x2_ave, label = "\$⟨x^2(t)⟩\$", color = "green")
    plot!(t, x4_ave, label = "\$⟨x^4(t)⟩\$", color = "red")

    pf = plot(p1, p2, p3, p4, layout = (2,2))
    savefig("StocEx1_moments.png")

    plot(t, σ2, ylabel = "\$ \\sigma ^2_{x^2}(t)\$", xlabel = "t", title="Variance of \$x^2\$",legend = false)
    savefig("StocEx1_sigmax2.png")
end

# Function for c)
function hist_plot(xt_data, times)
    interval = range(-2.25, 2.25, length = 50)
    plots = []
    for time in times
        t = time / 100
        p = histogram(xt_data[time,2:end], xlabel = "x($t)", ylabel = "Frequencies",
         bins = interval, normalize=:pdf, legend = false)
        
        plot!(p, xlabelfontsize=10, ylabelfontsize= 10, xticksfontsize=12, yticksfontsize= 12, size=(800, 600))
        push!(plots, p)
        
    end  
    plot(plots..., layout=(3,3))
    savefig("StocEx1_hist_time.png")
end

# Function for d)
# xt_data[:,2] --> First trayectory
function firstPassageTime(xt_data, xb, n)
    pass_t = Vector{Int64}()
    lim = size(xt_data)[1] # Dimension of each trayectory vector

    # i iterates over all the n= 5000 trayectories (+1 for the time)
    # j iterates over the points of each trayectory
    for i in 2:(n + 1)
        count = 0
        for j in 1:lim
            abs(xt_data[j,i]) < xb ? count+=1 : break
        end
        # Stores how many
        push!(pass_t, count)
    end

    return pass_t
end

# ----------------------------------------------------------------------------------------- #
# -----------------------          Ex1 a) "Trajectories:"           ----------------------- #
# ----------------------------------------------------------------------------------------- #

a, b, D = 4, 1, 0.01 # Parameters for the equation
x0 = 0.0 # Initial condition
Δt = 0.01 # Time between storing results
t_tot = 4 # Total time of the simulation
h = 0.001 # Time step
n = 5000 # Number of trayectories for integrating simo

# ⨰ = q(x,t) + g(x,t)ξ(t)
D = sqrt(D)
q(x, t) = a*x - b*x.^3 
g(x, t, D) = D
gprime(x, t) = 0

n_tot = t_tot / h # Total Number of iterations
n_stp = trunc(n_tot / (t_tot / Δt))   # Iterations between writing
n_wrt =  n_tot / n_stp  # Nº of times for repeating the writing loop

t = 0
x = fill(x0, n)
h_sqrt = sqrt(h)
file = open("trayectories_ex1.txt", "w")
for i in 1:n_wrt
    # Obtein the next point to write 
    x, t = milshtein(x, t, n_stp, h, h_sqrt, D, n)

    #Write results
    write(file, "$t")
    for num in x
        write(file, "\t$num")
    end
    write(file, "\n")
end
close(file)

# ----------------------------------------------------------------------------------------- #
# ----------------------- Ex1 b) "Transient anomalous fluctuations" ----------------------- #
# ----------------------------------------------------------------------------------------- #
xt_data = readdlm("trayectories_ex1.txt", '\t', Float64) # reading all the trayectories
moments_plot(xt_data)

# ----------------------------------------------------------------------------------------- #
# -----------------------  Ex1 c) "Probability density function:"   ----------------------- #
# ----------------------------------------------------------------------------------------- #s
times = Int64.(100 * [0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.5, 3, 4])
hist_plot(xt_data, times)

# ----------------------------------------------------------------------------------------- #
# ---------------------   Ex1 d) "First passage time distribution:"   --------------------- #
# ----------------------------------------------------------------------------------------- #
xb = 0.5
firstPassTime = firstPassageTime(xt_data, xb, n)
firstPassTime = firstPassTime * (h * n_stp) # Number of iterations to time

histogram(firstPassTime, xlabel = "Passage Time", ylabel = "Frequencies",
 normalize=:pdf, legend = false, title = "First Passage Time Distribution for \$ x_b = 0.5\$")

root = joinpath(@__DIR__, "images")
name = "\\StocEx1_PassTime.png"
savefig(root * name)

# ----------------------------------------------------------------------------------------- #
# ---------------------      Ex1 e) "Mean first passage time:"        --------------------- #
# ----------------------------------------------------------------------------------------- #
D_vect = [0.1, 0.01, 0.001, 0.0001]
D2_vect = sqrt.(D_vect)

time_ave = Vector{Float64}()
xb = 0.5

# We integrate 5000 trayectories for each value of D
for D_i in D2_vect
    count_vect = Vector{Float64}()
    for i in 1:5000
        x = [0]
        t = 0
        for j in 1:n_wrt
            x, t = milshtein(x, t, n_stp, h, h_sqrt, D_i, 1)
            if abs(x[1]) > xb
                break
            end
        end
        push!(count_vect, t)
    end
    push!(time_ave, mean(count_vect))
end

plot(log.(D_vect) , time_ave, marker = :circle, markersize = 3, 
 seriestype=:scatter, label = "\$ \\mu_t \$", xlabel = "Ln(D)", ylabel= "Mean First Passage Time (\$ \\mu_t \$)")
plot!(log.(D_vect), (-1/(2*a)) *log.(D_vect), label = "\$ -0.125lnD \$")

# Fit for the points
X = log.(D_vect)
Y = time_ave
cor(X,Y)
data = DataFrame(X=X, Y=Y)
modelo = lm(@formula(Y ~ X), data)
coeficientes = coef(modelo)
ftest(modelo.model)
X_pred = minimum(X .- 1):0.02:maximum(X .+ 1)
Y_pred = coeficientes[1] .+ coeficientes[2] * X_pred
plot!(X_pred, Y_pred, label="Linear Fit", color="lightgreen")
plot!(xlims=(-9.5,-1.5))

root = joinpath(@__DIR__, "images")
name = "\\StocEx1_Fit.png"
savefig(root * name)
end