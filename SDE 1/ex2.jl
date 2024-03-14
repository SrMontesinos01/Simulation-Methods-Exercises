module ex2

using DelimitedFiles, Plots, Measures, Statistics, ShiftedArrays, Trapz

# Function for integrating the Dif. Stoch. Eq.
function heun(x, y, t, h, h_sqrt, a, D, ϵ)
    uh_y = h_sqrt*randn()

    k_x = h.*q_x(x, y, ϵ)
    k_y = h.*q_y(x, a) .+ uh_y*g_y(D)    

    xk = x .+ k_x
    yk = y .+ k_y

    x = x .+ 0.5.*(k_x .+ h.*q_x(xk, yk, ϵ))
    y = y + 0.5.*(k_y .+ h.*q_y(xk, a) .+ uh_y*g_y(D))
    t = t + h
    return x, y, t
end

# Function for the Iterations for constant parameters
# If x0 and y0 are vectors, it will compute more than one trayectory
function simulation(x0, y0, h, h_sqrt, n_stp, n_wrt, a_s, D_s, ϵ_s)
    x = x0
    y = y0
    t = 0 
    file = open("trayectories_ex2.txt", "w")
    for i in 1:n_wrt
        #Write results
        write(file, "$t")
        for (num1, num2) in zip(x, y)
            write(file, "\t$num1\t$num2")
        end
        write(file, "\n")
    
        # Obtein the next point to write
        for j in 1:n_stp
            x, y, t = heun(x, y, t, h, h_sqrt, a_s, D_s, ϵ_s)
        end
    end
    close(file)
end

# Given a time series x(t), calculates the autocorrelation function
# The time should reach al least t = 80?
# tdif = h * n_stp
function AutroCorrFun(xt, s_max)
    xtnorm = xt .- mean(xt)
    xt_mean = mean(xtnorm.^2)

    # Since Δt = 0.1, multiply x10 to acces the vector element
    s_arr = [i for i in 0.1:0.1:s_max] * 10 
    Cs = []
    for s in s_arr
        # Prepare the shifted series
        xt_shift = ShiftedArrays.lead(xtnorm, Int64(trunc(s)))
        xt_shift = filter(!ismissing, xt_shift)
        xt_comp = xtnorm[1:length(xt_shift)]

        corr = mean(xt_comp .* xt_shift) / xt_mean
        push!(Cs, corr)
    end

    return s_arr/10, Cs
end

# Calculates the correlation time
function AutoCorrTime(corr, s) 
    corr = corr .^2
    return trapz(s, corr)
end

# Function for ploting various trayectories
function tray_plot(pos_data)
    n = size(pos_data)[2] - 1
    print(n)
    nums = [Int64(i) for i in 2:2:n]

    
    for num in nums
        k = num
        k2 = Int64(k/2)
        p_xy = plot(pos_data[:,k],pos_data[:,k+1], xlabel= "x(t)", ylabel = "y(t)", label="Trayectory")
        vline!([0], color=:black, linestyle=:solid, label=nothing)  # Vertical Line at x=0
        hline!([0], color=:black, linestyle=:solid, label=nothing)
        scatter!([x0[k2]], [y0[k2]], label="Initial Point", color = "red") # Plot the inital condition
        scatter!([x_st], [y_st], label="Fixed Point", color = "orange") # Plot the fixed point 

        p_xt = plot(pos_data[:,1],pos_data[:,k], xlabel= "t", ylabel = "x(t)", legend = false, label="x(t)")
        hline!([x_st], linestyle=:solid, label="Fixed Point", color = "orange")
        scatter!([0], [x0[k2]], label="Initial Point", color = "red")

        p_yt = plot(pos_data[:,1],pos_data[:,k+1], xlabel= "t", ylabel = "y(t)", legend = false, label="y(t)")
        hline!([y_st], linestyle=:solid, label="Fixed Point", color = "orange")
        scatter!([0], [y0[k2]], label="Initial Point", color = "red")

        # title = "\$D=$D\$" * " \$a=$a\$ " * " \$ ϵ=$ϵ\$ " * " \$ n_{stp} = $n_stp\$ " * " \$ n_{wrt} = $n_wrt\$ "
        p3 = plot(p_xy, p_xt, p_yt, layout =(1,3), size=(1800, 640), margin=15mm)
        plot(p3, legend =:outertop, legendfontsize=12, xtickfontsize=12, ytickfontsize=12, xlabelfontsize=12,
             ylabelfontsize=12, titlefontsize = 30)
    
        root = joinpath(@__DIR__, "images")
        name = "\\StocEx2_Tray"*string(num)*".png"
        savefig(root * name)
    end

end

# Function for plotting the Autocorrelation (diff. plots)
function corr_plot(s, c_arr, D_arr)
    p_arr = []
    for (corr, D) in zip(c_arr, D_arr)
        p = plot(s, corr, label="\$D = \$ $D" , xlabel = "Time Shift \$s\$", 
            ylabel = "Autocorrelation \$ C(s)\$")
        push!(p_arr, p)
    end


    pf = plot(p_arr..., layout = (2,3), ylims=(-0.80,1.2))

    root = joinpath(@__DIR__, "images")
    savefig(root * "\\StocEx2_corr_v2.png")
end

# ----------------- Funtions of the Dif. Eq. ----------------- # 
q_x(x, y, ϵ) = (x .- x .^3 ./ 3.0 .- y) / ϵ # Deterministic part of dx/dt
q_y(x, a) = a .+ x # Deterministic part of dy/dt
g_y(D) = D # Stochastic part of dy/dt

# ----------------------- Parameters ------------------------ # 
D, a, ϵ = [0.02, 0.04, 0.07, 0.1, 0.25, 0.9], 1.05, 0.01 # Parameters of the Dif. Eq.
h = 0.001 # Time step
h_sqrt, h2, h6 = sqrt(h), h/2.0, h/6.0

# t_tot = n_stp * n_wrt * h (Total time of the simulation)
n_stp = 100 # Iterations between writing
n_wrt = 5*500  # Nº of times for repeating the writing loop 

params = [D, a, ϵ, h, n_stp, n_wrt]
x_st = -a # x coord. for the fixed point (Deterministic part)
y_st = - a*(1 - a^2/3) # y coord. for the fixed point (Deterministic part)

# --------------------    Iterations    --------------------- # 
x0, y0 = [-2, 2, -1, 2], [-3, 3, -1.5, -3]
x0, y0 = [-1], [-1]
tau_arr = []
tau_std = []
corr_arr = []
reps = 1
s = 0
for D_i in D
    tau_reps = []
    for j in 1:reps
        simulation(x0, y0, h, h_sqrt, n_stp, n_wrt, a, D_i, ϵ)
        xy_data = readdlm("trayectories_ex2.txt", '\t', Float64)
        yt = xy_data[:,3]

        s, corr = AutroCorrFun(yt, 40)
        push!(corr_arr, corr)
        tau = AutoCorrTime(corr, s)
        push!(tau_reps, tau)
    end 

    push!(tau_arr, mean(tau_reps))
    push!(tau_std, std(tau_reps))
end

xy_data = readdlm("trayectories_ex2.txt", '\t', Float64)
tray_plot(xy_data) # For storing all the trayectories in xy_data
corr_plot(s, corr_arr, D)
plot(D, tau_arr, yerror = tau_std, xlabel= "D", ylabel = "Autocorrelation Time",
     legend=:none, marker=:circle)

root = joinpath(@__DIR__, "images")
name = "\\StocEx2_corrTime.png"
savefig(root * name)

D
tau_arr
# -------------- Testing Stuff -------------- #
simulation(x0, y0, h, h_sqrt, n_stp, n_wrt, a, 0.02, ϵ)
xy_data = readdlm("trayectories_ex2.txt", '\t', Float64)
yt = xy_data[:,3]
plot(yt)
s, corr = AutroCorrFun(yt, 40)

xtnorm = yt .- mean(yt)
xt_mean = mean(yt .*yt)
xt_shift = ShiftedArrays.lead(xtnorm, Int64(trunc(70)))
xt_shift = filter(!ismissing, xt_shift)
xt_comp = xtnorm[1:length(xt_shift)]
plot(xtnorm[1:300])
plot!(xt_shift[1:300])
plot!(xt_comp[1:300])
corr = mean(xt_comp .* xt_shift) / xt_mean

plot(s, corr, xlabel = "Time Shift \$s\$", ylabel = "Autocorrelation \$ C(s)\$")
plot(xy_data[:,2],xy_data[:,3])
plot(xy_data[:,1],xy_data[:,3]) 
# -------------------------------------------- #



# Plotting all the trayectories at once
plot(xy_data[:,2],xy_data[:,3])
plot!(xy_data[:,4],xy_data[:,5])
plot!(xy_data[:,6],xy_data[:,7])
plot!(xy_data[:,8],xy_data[:,9])
scatter!([x0], [y0], label="Initial Point") # Plot the inital condition
scatter!([x_st], [y_st], label="Initial Point") # Plot the fixed point
vline!([0], color=:black, linestyle=:solid, label=nothing)  # Línea vertical en x=0
hline!([0], color=:black, linestyle=:solid, label=nothing)  # Línea horizontal en y=0

end
