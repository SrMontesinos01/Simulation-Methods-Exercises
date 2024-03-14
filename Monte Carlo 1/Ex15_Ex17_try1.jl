# ------------------------------------------------------------------- #
# Hit and Miss Monte Carlo method
# "Experiment": Calculating the integral using M random uniform points
# ------------------------------------------------------------------- #
using Distributions, Plots, Statistics, Random
using GLM, DataFrames

# Function that we want to integrate
fun(x::Real) = sqrt(1 - x^2)

# Function for the Hit and Miss MC method
function hit_and_miss(a , b, c, M; method = 1)
    cont = 0
    for i in 1:M
        if method == 1
            x = rand(Uniform(a, b))
            y = rand(Uniform(0, c))
        elseif method ==2
            x = rand(MersenneTwister(), Uniform(a, b))
            y = rand(MersenneTwister(), Uniform(0, c))
        else
            x = rand(RandomDevice(), Uniform(a, b))
            y = rand(RandomDevice(), Uniform(0, c))
        end

        if fun(x) > y 
            cont += 1 
        end
    end

    prob = cont / M
    int = (b - a) * c * prob
    sigma = (b - a)c * sqrt( (1 - prob)prob / M)

    return int, sigma
end

# Macro for calculating the time
@time begin
    
    a , b, c = 0, 1, 1  # Limits for the surface to explore
    # Set of Numbers of points for the method
    M_list = [10^x for x in 1:6] # Ex15
    # M_list = collect(range(10, 10^4, 15)) # Ex15
    M_list = [100] # Ex17

    μ_val_list = [] # For storing the mean of real erros
    error_list = []
    μ_σ_list = [] # For storing the mean of sigma 
    σ_error_list = []
    percentage = [0,0,0,0]
    reps = 10^4 # number of times that we repeat each "experiment"

    for M in M_list
        val_list = [] # For calculating the mean real error from experiments
        var_list = [] # For calculating the mean sigma from experiments
        for i in 1:reps
            int, sigma = hit_and_miss(a , b, c, M; method = 3)

            # Calculate the real error and append it 
            sol = π/4 # Exact Solution
            real_error = abs(int - sol)
            push!(val_list, real_error)
            push!(var_list, sigma)

            if real_error < sigma
                percentage[1] +=1 
            elseif real_error < 2 * sigma
                percentage[2] +=1 
            elseif real_error < 3 * sigma
                percentage[3] +=1
            elseif real_error > 4 * sigma
                percentage[4] +=1
            end
            
        end 
        
        percentage = (percentage / reps) * 100
        percentage[2] += percentage[1]
        percentage[3] += percentage[2]

        # Calculate a mean real error and its standard deviation
        push!(μ_val_list, mean(val_list))
        push!(error_list, std(val_list))

        push!(μ_σ_list, mean(var_list))
        push!(σ_error_list, std(var_list))
    end

    print(percentage)
    sum(percentage)
end

# ------------------------------------------------------------------- #
# --------------------- Ploting results sigma-M --------------------- #
# ------------------------------------------------------------------- #
plot1 = plot(M_list, μ_σ_list,
            yerror = σ_error_list ,
            title = " σ vs M" ,
            xscale=:log10, 
            msc=:auto, 
            xlabel = "M (number of points for MC method)", 
            ylabel = "σ",
            label = "σ",
            marker = :circle,
            markersize = 3,
            seriestype=:scatter)

m_σ(x) = 0.41 / sqrt(x)
plot1 = plot!(m_σ,
            xscale=:log10,
            label = "0.41/√(M)")

# For the fit
X = 1 ./ sqrt.(M_list)
Y =  1.0 .* μ_σ_list # Must be Float64
cor(X,Y)
data = DataFrame(X=X, Y=Y)
modelo = lm(@formula(Y ~ X), data)
coeficientes = coef(modelo)
ftest(modelo.model)
X_pred = minimum(X):0.02:maximum(X .+ 0.02)
Y_pred = coeficientes[1] .+ coeficientes[2] * X_pred

plot2 =plot(X, Y, 
            yerror = σ_error_list ,
            title = " σ vs 1 / √(M)" , 
            msc=:auto, 
            xlabel = "1/√(M)", 
            ylabel = "σ",
            label = "σ",
            marker = :circle,
            markersize = 3,
            seriestype=:scatter
 )

plot2 = plot!(X_pred, Y_pred, label="Linear Fit", color="lightgreen")
plotd = plot(plot1, plot2, layout=(1, 2))
savefig(plotd,"sigma_fit.png")

# ------------------------------------------------------------------- #
# ------------------- Ploting results real error - M ---------------- #
# ------------------------------------------------------------------- #
plot_error = plot(M_list, μ_σ_list,
 yerror = σ_error_list ,
 title = "" ,
 xscale=:log10, yscale=:log10,
 msc=:auto, 
 label = "σ",
 marker = :circle,
 markersize = 3,
 seriestype=:line,
 )

plot_error = plot!(M_list, μ_val_list,
 yerror = error_list ,
 xscale=:log10, yscale=:log10,
 msc=:auto, 
 xlabel = "M (number of points for MC method)", 
 ylabel = "Error",
 label = "Real Error",
 marker = :circle,
 markersize = 3,
 seriestype=:line,
 color = "darkorange")



X = 1.0 * log10.(M_list)
Y = 1.0 * log10.(μ_val_list) # Esto pq tiene que ser Float64
cor(X,Y)
plot(X, Y, label="Linear Fit", color="lightgreen")


data = DataFrame(X=X, Y=Y)
modelo = lm(@formula(Y ~ X), data)
coeficientes = coef(modelo)
ftest(modelo.model)
X_pred = minimum(X):0.02:maximum(X .+ 0.02)
Y_pred = coeficientes[1] .+ coeficientes[2] * X_pred

plot_fit2 = plot(X_pred, Y_pred,
    xlabel = "log_10(M) ", 
    ylabel = "log_10(Real Error)",
    label="Linear Fit",
    color="lightgreen")

plot_fit2 = plot!(X, Y,
    xlabel = "log_10(M) ",
    msc=:auto,  
    ylabel = "log_10(Real Error)",
    label="Real Error",
    marker = :circle,
    markersize = 3,
    seriestype=:line,
    color= "darkorange")    
    

plotd = plot(plot_error, plot_fit2, layout=(1, 2))
savefig(plotd,"real_error_fit.png")
