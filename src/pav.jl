## pav(y) returns the isotonic regression of the predictor y.
## it is slow, and possibly wrong.   
function pav{T<:Real}(y::Vector{T})
    yy = convert(Vector{Float64}, y)
    n = length(y)
    i::Int=1
    while i<n
        while yy[i+1] >= yy[i] 
            if i<n-1
                i += 1
            else
                return yy
            end
#            println("i ", i)
        end
        ie = i+1
#	println("range ", i, " ", ie)
        yy[i:ie] = mean(y[i:ie])
#       println("yy+ ", yy')
        while i>1 && yy[i-1] > yy[i]
            i -= 1
            yy[i:ie] = mean(y[i:ie])
        end
#        println("yy- ", yy')
    end
    yy
end
