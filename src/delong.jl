using Distributions

function ψ(x,y)
    return if x > y 
        1
        elseif x < y
        0
        else
        1/2
    end
end
    
function V₁₀(x::Float64; Y)
    sum(ψ.(x,Y)) / length(Y)
end

function V₀₁(y::Float64; X)
    sum(ψ.(X, y)) / length(X)
end
    
function V₁₀(X::Vector{N},Y::Vector{N}) where {N<:Real}
    V₁₀.(X, Y=Y)
end

function V₀₁(X::Vector{N},Y::Vector{N}) where {N<:Real}
    V₀₁.(Y, X=X)
end
    
function S₁₀(tars, nontars, θs)
    Xs = hcat(map(zip(tars, nontars, θs)) do (tar, nontar, θ)
        V₁₀.(tar, Y=nontar) .- θ
    end...)
    cov(Xs)
end
    
function S₀₁(tars, nontars, θs)
    Xs = hcat(map(zip(tars, nontars, θs)) do (tar, nontar, θ)
        V₀₁.(nontar, X=tar) .- θ
    end...)
    cov(Xs)
end
        
function delong_test(tar1, nontar1, tar2, nontar2)
    

    # Adjust the lengths in case of different lengths
    m = minimum(length.([tar1, tar2]))
    n = minimum(length.([nontar1, nontar2]))
        
    tars = [tar1[1:m], tar2[1:m]]
    nontars = [nontar1[1:n], nontar2[1:n]]
        
    θ₁ = auc(tar1, nontar1)
    θ₂ = auc(tar2, nontar2)
    θ  = [θ₁, θ₂]
        
    l = [1, -1]
    
    S₁₀(tars, nontars, θ)
    S₀₁(tars, nontars, θ)
    Delon_stat = abs((θ₁ - θ₂)) /  √(l' * (S₁₀(tars, nontars, θ)/m + S₀₁(tars, nontars, θ)/n) * l)
    
    1 - cdf(Normal(0,1), Delon_stat) + cdf(Normal(0,1), -Delon_stat)
end