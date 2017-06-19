



function discretiseParametrical(h, s, t)
    H = [h(σ,τ) for τ in t, σ in s]

    X = [ C[1] for C in H ]
    Z = [ C[2] for C in H ]

    TwoVector(X, Z)
end

function topContainer(s,t)
    X = zeros(length(t), length(s))
    Z = zeros(length(t), length(s))

    X[1,:] = s

    TwoVector(X, Z)
end

function phiContainer(s,t)
    zeros(length(t), length(s))
end



function bottomNormal(C::TwoVector)
    curveNormal(C)
end

function surfaceNormal(C::TwoVector)
    -1. * curveNormal(C)
end



#CONVERSIONS

function sToα(C::TwoVector)
    s = size(C.X)

    X = zeros(s[1],s[2]+1)
    Z = zeros(s[1],s[2]+1)

    for i in 1:s[1],j in 2:s[2]
        X[i,j] = ( C.X[i,j] + C.X[i,j-1] )/2
        Z[i,j] = ( C.Z[i,j] + C.Z[i,j-1] )/2

    end

    TwoVector(X,Z)

end

function αTos(C::TwoVector)
    s = size(C.X)

    X = zeros(s[1],s[2]-1)
    Z = zeros(s[1],s[2]-1)

    for i in 1:(s[1]-1),j in 1:(s[2]-1)
        X[i,j] = ( C.X[i,j+1] + C.X[i,j] )/2
        Z[i,j] = ( C.Z[i,j+1] + C.Z[i,j] )/2
    end

    TwoVector(X,Z)
end

function sToαOfLine(L)
    l = length(L)
    La = Array(typeof(L[1]), l+1)
    La[1] = zeros(length(L[1]))'
    La[end] = zeros(length(L[1]))'

    for j in 2:l
        La[j] = ( L[j] + L[j-1] )/2
    end

    La
end

function αTosOfFloatLine(L)
    l = length(L)
    La = zeros(l-1)

    for j in 1:l-1
        La[j] = ( L[j] + L[j+1] )/2
    end

    La
end




function solveLinearSystem(ϕ1, ∇nϕ2, C1,C2,N1,N2)

    A11 = A(1, 1, C1, C2, N1, N2)
    A12 = A(1, 2, C1, C2, N1, N2)
    A21 = A(2, 1, C1, C2, N1, N2)
    A22 = A(2, 2, C1, C2, N1, N2)

    B11 = B(1, 1, C1, C2, N1, N2)
    B12 = B(1, 2, C1, C2, N1, N2)
    B21 = B(2, 1, C1, C2, N1, N2)
    B22 = B(2, 2, C1, C2, N1, N2)

    L = [ -B11 A12; -B21 A22-0.5*eye(ns,ns)]
    R = [ 0.5*eye(ns,ns)-A11 B12; -A21 B22 ]

    #System: Lu = Ru0

    u0 = vec([ϕ1 ∇nϕ2])

    Ru0 = R*u0

    u = L\Ru0

    ( u[1:ns], u[(ns+1):end] )
end





function gradOfLine(L)
    l = length(L)

    dL = Array(typeof(L[1]), l)

    for j in 1:l
        if j == 1
            dL[j] = ( -3*L[j] + 4*L[j+1] - L[j+2] )/(2*ds)
        elseif j == l
            dL[j] = ( L[j-2] - 4*L[j-1] + 3*L[j] )/(2*ds)
        else
            dL[j] = ( L[j+1] - L[j-1] )/(2*ds)
        end
    end

    dL
end

function tangentOfLine(L)
    D = diff(L)
    e = [norm(d) for d in D]

    D./e
end

function surfaceNormalOfLine(L)
    D = diff(L)

    N = [ [C[2] -C[1]] for C in D]
    e = [norm(d) for d in D]

    -1*N./e
end
