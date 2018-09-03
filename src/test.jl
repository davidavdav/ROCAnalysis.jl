t1 = ([1], [0])
t2 = ([0], [1])
t3 = ([0], [-1, 1])
t4 = ([-1, 1], [0])

t5 = ([1:2], [0:1])

## pav(y) returns the isotonic regression of the predictor y.
## after W.J. Wilbur, L. Yeganova and W. Kim (2005)
function pav(y::Vector{T}) where T<:Real
    yy = convert(Vector{Float64}, y)
    N = length(y)
    F = [1:N]+1                 # forward pointer to intervals
    B = F-2                     # backward pointer to intervals
    Q1 = [1:N]                # queues
    Q2 = Int[]
    flag = false
    while length(Q1)>0
        m = 1
        while length(Q1)>0
            k = shift!(Q1)
            if m <= k
                flag=false
                println(k, " ", F[k])
                while F[k] != N && yy[k]> yy[F[k]]
                    ## pool intervals beginning at k and F[k]
                    first = k>0 ? F[k-1] : 1
                    pooled = vcat(first:F[k], F[F[k]-1]:F[F[k]])
                    println("pool ", pooled')
                    yy[pooled] = mean(y[pooled])
                    u = F[F[k]]
                    F[k] = u
                    if u<=N B[u] = k end
                    m = u
                    flag = true
                end
            end
            if flag && B[k]>=0
                push!(Q2, B[k])
            end
        end
        h = Q2; Q2 = Q1; Q1 = h
    end
    yy
end
