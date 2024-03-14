using Plots, GLM, DataFrames
# Datos de ejemplo
X = [1.0, 2.0, 3.0, 4.0, 5.0]
Y = [2.0, 3.0, 4.0, 4.5, 6.0]

# Crear un DataFrame con los datos
data = DataFrame(X=X, Y=Y)

# Ajuste lineal
modelo = lm(@formula(Y ~ X), data)

# Obtener coeficientes del modelo
coeficientes = coef(modelo)

# Crear un rango de valores para X
X_pred = minimum(X):0.1:maximum(X)

# Calcular los valores predichos usando el modelo ajustado
Y_pred = coeficientes[1] .+ coeficientes[2] * X_pred

# Crear el gráfico
using Plots
scatter(X, Y, label="Datos", legend=:topright)
plot!(X_pred, Y_pred, label="Ajuste Lineal", xlabel="X", ylabel="Y")

