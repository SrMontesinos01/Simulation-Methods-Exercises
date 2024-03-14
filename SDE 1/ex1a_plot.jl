using DelimitedFiles

# ----------------------------------------------------------------------------------------- #
# -----------------------           Ex1 a) "Trajectories"           ----------------------- #
# ----------------------------------------------------------------------------------------- #
function ploting(matrix)
    lim = size(matrix)[2]
    t = matrix[:,1]
    x = matrix[:,2]
    plot(t, x, legend = false, 
        xlabel = "t", ylabel = "x(t)", title = "Milshtein Algorithm: Integrated Trajectories")

    for i in 3:lim
        x = matrix[:,i]
        plot!(t, x)
        gui()
    end

    return current()
end

# Read the wholw file 
matrix = readdlm("trayectories_ex1.txt", '\t', Float64)

# Make the plot
p = ploting(matrix)
savefig("StocEx1_20tray.png")

