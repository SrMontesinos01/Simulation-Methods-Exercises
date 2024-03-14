import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.widgets import Button, Slider

# Define el sistema de ecuaciones
def system(t, u, e, a):
    
    x, y = u
    dxdt = (x - (x**3) /3.0 - y) / e
    dydt = x + a
    return [dxdt, dydt]

# Define una función para calcular las nullclinas
def nullx(x, mu): 
    value = mu * (x**3 / 3.0 - x)
    return value
    
e = 0.01 # Parametro del sistema
a = 1.05  # Parametro del sistema

u0 = [3, -10] # Condiciones iniciales
t_span = (0, 150) # Rango de tiempo para la integración
m = 1 # Margen para la trayectoria al hacer el plot

# Integra para obtener la trayectoria (x,y)
sol = solve_ivp(system, t_span, u0,
                t_eval=np.linspace(t_span[0], t_span[1], 1000),
                args=(e, a))
x_values = sol.y[0]
y_values = sol.y[1]

# Crea una malla de valores para el campo vectorial
x_grid, y_grid = np.meshgrid(np.linspace(min(x_values) - m, max(x_values) + m, 20 + 4*m), 
                             np.linspace(min(y_values) - m, max(y_values) + m, 20 + 4*m))

# Calcula las derivadas dx/dt y dy/dt en la malla
dxdt_grid, dydt_grid = system(0, [x_grid, y_grid], e, a)

# Calcula la magnitud del campo vectorial
magnitude = np.sqrt(dxdt_grid**2 + dydt_grid**2)


# Escala logarítmica para reducir las flechas más largas
scale_factor = 5  # Factor de escala
U = dxdt_grid
V = dydt_grid
scaled_U = np.log(1 + np.abs(U)) * np.sign(U) / scale_factor
scaled_V = np.log(1 + np.abs(V)) * np.sign(V) / scale_factor

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++++++++++++++++++++++++ PLOT STUFF +++++++++++++++++++++++++ #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Grafica el campo vectorial con un mapa de color
plt.figure(figsize=(8, 6))
plt.quiver(x_grid, y_grid, scaled_U, scaled_V, magnitude, cmap='hot', scale=20, alpha=0.80)
plt.colorbar(label='Magnitud del campo vectorial')

# Grafica las nulclinas
# x_null = np.linspace(-5, 5, 100) 
# y_null = nullx(x_null, mu)
# plt.plot(x_null, y_null, linestyle = "dashed", linewidth=2)
# plt.axvline(a, linestyle = "dashed", color='red', linewidth=2) # Añade una línea vertical en x=0 (eje y)

# Trazar la línea de puntos como cruces verdes
plt.plot(x_values, y_values, 'g-', linewidth=2, label='Solución')
plt.scatter(u0[0], u0[1], color='green', marker='o', label='x0') # Cond. Inicial

plt.xlim(min(x_values) - m, max(x_values) + m)  # Límites del eje x
plt.ylim(min(y_values) - m, max(y_values) + m)  # Límites del eje y

plt.axhline(0, color='black', linewidth=2) # Añade una línea horizontal en y=0 (eje x)
plt.axvline(0, color='black', linewidth=2) # Añade una línea vertical en x=0 (eje y)


plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.show()

