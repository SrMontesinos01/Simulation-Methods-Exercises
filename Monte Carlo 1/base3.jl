function decimal_to_base3(number, n)
    # Inicializa un vector de dimensión n lleno de ceros
    base3_vector = zeros(Int, n)
    
    # Realiza la conversión a base 3 y almacena los dígitos en el vector
    i = 1
    while number > 0 && i <= n
        remainder = number % 3
        base3_vector[i] = remainder
        number = div(number, 3)
        i += 1
    end
    
    return reverse(base3_vector) .- 1
end

# Función para quitarle 1 a cada elemento del vector
function subtract_one_from_base3(base3_vector)
    return base3_vector .- 1
end

# Ejemplo de uso
decimal_number = 17
dimension_n = 6

base3_result = decimal_to_base3(decimal_number, dimension_n)
println("Número $decimal_number en base 3 con dimensión $dimension_n: ", base3_result)

result_after_subtraction = subtract_one_from_base3(base3_result)
println("Resultado después de quitarle 1 a cada elemento: ", result_after_subtraction)
