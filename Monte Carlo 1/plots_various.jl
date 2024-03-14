using Plots

# Sample data
x = 1:10
y1 = sin.(x)
y2 = cos.(x)
y3 = tan.(x)

# Create a 2x2 grid of subplots
plot1 = plot(x, y1, label="sin(x)", xlabel="x", ylabel="y")
plot1 = plot!(x, y2, label="a(x)", xlabel="x", ylabel="y")
plot2 = plot(x, y2, label="cos(x)", xlabel="x", ylabel="y")
plot3 = plot(x, y3, label="tan(x)", xlabel="x", ylabel="y")

# Combine the subplots into a single plot
plot(plot1, plot2, plot3, layout=(1, 3))
