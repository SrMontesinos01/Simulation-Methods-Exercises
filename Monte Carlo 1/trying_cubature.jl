using Cubature

function fun(x::Vector{Float64})
    s = 0.0  
    
    for i in eachindex(x)
        
        i == lastindex(x) ? s += x[1] * x[i] : s += x[i] * x[i + 1]
    end

    value = exp(-sum(x .^ 2)) * cos(s)^2
    return value
end

a = -1
b = 1

inf = [a for i in 1:5]
sup = [b for i in 1:5]
intm, errm = hcubature(fun, inf, sup)

f(x) = 2*x
int, err = hquadrature(f, 0.0, 1.0)