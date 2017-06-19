
import Base
import Base.LinAlg: norm, *, ./, transpose
#'
type TwoVector
    X::Matrix{Float64}
    Z::Matrix{Float64}
end


function *(a::Number, C::TwoVector)
    TwoVector(a*C.X,a*C.Z)
end
function ./(C::TwoVector, D::Array)
    TwoVector(C.X./D, C.Z./D)
end

function norm(V::TwoVector)
   sqrt(V.X.^2 + V.Z.^2)
end

function transpose(C::TwoVector)
   TwoVector(C.X', C.Z')
end






function spacediff(C::TwoVector)
    dX = diff(C.X,2)
    dZ = diff(C.Z,2)

    TwoVector(dX,dZ)
end

function timediff(C::TwoVector)
    dX = diff(C.X,1)
    dZ = diff(C.Z,1)

    TwoVector(dX,dZ)
end



#Difference between adjacent points on curve divided by parametrisation interval
function grad(C::TwoVector)
    s = size(C.X)

    dX = zeros(s)
    dZ = zeros(s)

    for i in 1:s[1],j in 1:s[2]
        if j == 1
            dX[i,j] = ( -3*C.X[i,j] + 4*C.X[i,j+1] - C.X[i,j+2] )/(2*ds)
            dZ[i,j] = ( -3*C.Z[i,j] + 4*C.Z[i,j+1] - C.Z[i,j+2] )/(2*ds)
        elseif j == s[2]
            dX[i,j] = ( C.X[i,j-2] - 4*C.X[i,j-1] + 3*C.X[i,j] )/(2*ds)
            dZ[i,j] = ( C.Z[i,j-2] - 4*C.Z[i,j-1] + 3*C.Z[i,j] )/(2*ds)
        else
            dX[i,j] = ( C.X[i,j+1] - C.X[i,j-1] )/(2*ds)
            dZ[i,j] = ( C.Z[i,j+1] - C.Z[i,j-1] )/(2*ds)
        end
    end

    TwoVector(dX,dZ)
end
function timeDerivative(C::TwoVector)
    s = size(C.X)

    dX = zeros(s)
    dZ = zeros(s)

    for i in 1:s[1],j in 1:s[2]
        if i == 1
            dX[i,j] = ( -3*C.X[i,j] + 4*C.X[i+1,j] - C.X[i+2,j] )/(2*dt)
            dZ[i,j] = ( -3*C.Z[i,j] + 4*C.Z[i+1,j] - C.Z[i+2,j] )/(2*dt)
        elseif i == s[1]
            dX[i,j] = ( C.X[i-2,j] - 4*C.X[i-1,j] + 3*C.X[i,j] )/(2*dt)
            dZ[i,j] = ( C.Z[i-2,j] - 4*C.Z[i-1,j] + 3*C.Z[i,j] )/(2*dt)
        else
            dX[i,j] = ( C.X[i+1,j] - C.X[i-1,j] )/(2*dt)
            dZ[i,j] = ( C.Z[i+1,j] - C.Z[i-1,j] )/(2*dt)
        end
    end

    TwoVector(dX,dZ)
end


function curveTangent(C::TwoVector)
    D = spacediff(C)
    e = norm(D)

    D./e
end

#TODO: Implement signs of paths, DONE: in specific implementations for the different curves
function curveNormal(C::TwoVector)
    D = spacediff(C)
    N = TwoVector(D.Z,-D.X)
    e = norm(D)

    N./e
end
