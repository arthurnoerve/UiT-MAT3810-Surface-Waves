

#Matrix elements



#Returns NxN matrices for the A and B coefficient matrices
function A(i::Int, j::Int, C1, C2, N1, N2)
    e1 = [ norm(p) for p in diff(C1) ]
    e2 = [ norm(p) for p in diff(C2) ]

    if i == 1 && j == 1
        function element11(k,l)
            if i==j || k==l
                return 0
            else
                v = C1[l] - C1[k]
                n = N1[l]
                return e1[l]/(2π)*( (dot(v,n))/norm(v)^2 )
            end
        end
        return [element11(k,l) for k in 1:ns, l in 1:ns]
    elseif i == 1 && j == 2
        function element12(k,l)
            if i==j || k==l
                return 0
            else
                v = C2[l] - C1[k]
                n = N2[l]
                return e2[l]/(2π)*( (dot(v,n))/norm(v)^2 )
            end
        end
        return [element12(k,l) for k in 1:ns, l in 1:ns]
    elseif i == 2 && j == 1
        function element21(k,l)
            if i==j || k==l
                return 0
            else
                v = C1[l] - C2[k]
                n = N1[l]
                return e1[l]/(2π)*( (dot(v,n))/norm(v)^2 )
            end
        end
        return [element21(k,l) for k in 1:ns, l in 1:ns]
    elseif i == 2 && j == 2
        function element22(k,l)
            if i==j || k==l
                return 0
            else
                v = C2[l] - C2[k]
                n = N2[l]
                return e2[l]/(2π)*( (dot(v,n))/norm(v)^2 )
            end
        end
        return [element22(k,l) for k in 1:ns, l in 1:ns]
    end



end





function B(i::Int, j::Int, C1, C2, N1, N2)
    e1 = [norm(p) for p in diff(C1) ]
    e2 = [norm(p) for p in diff(C2) ]

    if i == 1 && j == 1
        function element11(k,l)
            if i==j || k==l
                return e1[l]/(2π) * log(e1[l]) - 1
            else
                v = C1[l] - C1[k]
                return e1[l]/(2π) * log( norm(v) )
            end

        end
        return [element11(k,l) for k in 1:ns, l in 1:ns]
    elseif i == 1 && j == 2
        function element12(k,l)
            if i==j || k==l
                return e2[l]/(2π) * log(e2[l]) - 1
            else
                v = C2[l] - C1[k]
                return e2[l]/(2π) * log( norm(v) )
            end

        end
        return [element12(k,l) for k in 1:ns, l in 1:ns]
    elseif i == 2 && j == 1
        function element21(k,l)
            if i==j || k==l
                return e1[l]/(2π) * log(e1[l]) - 1
            else
                v = C1[l] - C2[k]
                return e1[l]/(2π) * log( norm(v) )
            end

        end
        return [element21(k,l) for k in 1:ns, l in 1:ns]
    elseif i == 2 && j == 2
        function element22(k,l)
            if i==j || k==l
                return e2[l]/(2π) * log(e2[l]) - 1
            else
                v = C2[l] - C2[k]
                return e2[l]/(2π) * log( norm(v) )
            end

        end
        return [element22(k,l) for k in 1:ns, l in 1:ns]
    end

end
