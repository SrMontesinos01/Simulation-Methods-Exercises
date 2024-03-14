using Distributions, Plots, Statistics, Random

function fun(x::Vector{Float64})
    s = 0.0  
    for i in eachindex(x)
        i == lastindex(x) ? s += x[1] * x[i] : s += x[i] * x[i + 1]
    end

    value = exp(-sum(x .^ 2)) * cos(s)^2
    return value
end

# Function for the Hit and Miss MC method
function hit_and_miss(a , b, c, M, n; method = 1)
    cont = 0
    for i in 1:M
        if method == 1
            x = [rand(Uniform(a, b)) for i in 1:n]
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
    int = (b - a)^n * c * prob
    sigma = (b - a)^n * sqrt( (1 - prob)prob / M) * c

    return int, sigma
end

#Function for the Uniform Sampling
function uniform_sampling(a, b, M, n)
    G_list = []
 
    for i in 1:M
        x = [rand(Uniform(a, b)) for i in 1:n]
        push!(G_list, fun(x) )
    end

    G_list = (b -a)^n .* G_list
    int = sum(G_list) / M
    sigma = sqrt( sum(G_list.^2) / M - (sum(G_list) / M)^2 ) / sqrt(M)
    
    return int, sigma
end

# Auxiliar function for Simpson implementation
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
    
    return reverse(base3_vector).- 1
end

# Function for N-dim Simpson Rule
function ND_simp(a, b, n)

    limit = 3^n -1
    suma = 0
    for i in 0:(limit)
        n_vect = decimal_to_base3(i, n)
        expo= n - sum(n_vect.^2)
        suma +=  fun(1.0 * n_vect) * 4^(expo)
    end

    suma = suma * (1/3)^n
    return  suma 
end

a , b, c = -1, 1, 1  # Limits for the surface to explore
M , n = 10^6, 10 # Number of point and dimesion of the problem

@time begin
    hit_and_miss(a, b, c, M, n)
end

@time begin
    uniform_sampling(a, b, M, n)
end

@time begin
    ND_simp(a, b, n)
end
