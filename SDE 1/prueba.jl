using Plots
s1 = scatter(1:10, 1:10, legend = false)
s2 = scatter(1:10, 1:10, legend = false)
s3 = scatter(1:10, 1:10, legend = false)
s4 = scatter(1:10, 1:10, legend = false)
legend = plot([0 0 0 0], showaxis = false, grid = false, label = ["s1" "s2" "s3" "s4"])
plot(s1, s2, s3, s4, legend, layout = @layout([[A B; C D] E{.1w}]))

using MakieLayout
using ColorSchemes

using AbstractPlotting
using ShiftedArrays
# Función para realizar un desplazamiento circular en un vector
function circular_shift(vector, s)
    n = length(vector)
    s = s % n  # Asegurarse de que s esté en el rango [0, n-1]
    return [vector[(i - s) % n + 1] for i in 1:n]
end


# Supongamos que tenemos una serie temporal y(t) almacenada en un vector
y = [1, 2, 3, 4, 5, 6, 7, 8, 9]
y = ShiftedArray(y,2)
typeof(y[1])
y = ShiftedArrays.lead(y, 2)
y = ShiftedArrays.circshift(y,2)
filtered_vec = filter(!ismissing, y)
lead()
# Especifica el desplazamiento s que deseas aplicar
s = 2  # Cambia el valor de s según tu necesidad

# Calcula la serie temporal desplazada y(t + s) utilizando la función personalizada
y_desplazada = circular_shift(y, s)

# y_desplazada contendrá la serie temporal desplazada y(t + s)
println(y_desplazada)


# Borrado del ex2
s_arr = [i for i in 0:0.1:40]*100
xt = [i for i in 1:8000]
n = 300
xt_shift = ShiftedArrays.lead(xt, Int64(700))
xt_shift = filter(!ismissing, xt_shift)
xt_shift = xt_shift[1:4000]