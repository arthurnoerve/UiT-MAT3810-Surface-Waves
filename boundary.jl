

function boundaryIntegral(L,ns,T,nt,h)
    ds = 2L/ns
    dt = T/nt

    #Grid discretisation
    s = collect(linspace(-L,L,ns+1))
    smid = [-L + ds*(k-1/2) for k in 1:ns]

    t = collect(linspace(0,T,nt))

    #CONTAINERS

    #bottom
    ϕ2 = phiContainer(smid,t)  #6.10
    C2 = discretiseParametrical(h,s,t)
    n2 = bottomNormal(C2) #Normal to bottom
    p2 = curveTangent(C2) #Tangent to bottom


    #surface
    ϕ1 = phiContainer(smid,t)  #6.10
    C1 = topContainer(s,t)
    n1 =  TwoVector(zeros(ns,ns),zeros(ns,ns))#Normal to bottom
    p1 =  TwoVector(zeros(ns,ns),zeros(ns,ns))#Tangent to bottom




    #PRECOMPUTATION
    ∇ϕ2 = αTos(timeDerivative(C2)) #6.5
    ∇n2 = grad(n2)

    #Combine ∇ϕ2 and ∇n2 from above to form ∇nϕ2
    ∇nϕ2 =  ∇n2.X .* ∇ϕ2.X + ∇n2.Z .* ∇ϕ2.Z #7.8



    #COMPUTATION

    #EULER LOOP
    S0 = [ C1.X[1,:], C1.Z[1,:], ϕ1[1,:] ]
    S = Array(Array{Array{Float64,1},1},nt)
    S[1] = S0

    for k in 1:(length(t)-1)
        #Extract
        ϕ1E = S[k][3]
        ∇nϕ2E = ∇nϕ2[k,:]
        C1E = [ [S[k][1][i] S[k][2][i]] for i in 1:ns+1 ]
        C2E = [ [C2.X[k,i] C2.Z[k,i]] for i in 1:ns+1 ]
        n2E = [ [n2.X[k,i] n2.Z[k,i]] for i in 1:ns ]

        #Compute normal and tangent
        n1E = surfaceNormalOfLine(C1E)
        p1E = tangentOfLine(C1E)


        (∇nϕ1E, ϕ2E) = solveLinearSystem(ϕ1E,∇nϕ2E,C1E,C2E,n1E,n2E) #7.3

        ∇pϕ1E = gradOfLine(ϕ1E) #6.8

        ∇n1E = gradOfLine(n1E)
        ∇p1E = gradOfLine(p1E)

        ∇ϕ1E_mid = ∇nϕ1E .* ∇n1E + ∇pϕ1E .* ∇p1E #7.7
        ∇ϕ1E = sToαOfLine(∇ϕ1E_mid)


        ϕX_mid = zeros(length(∇ϕ1E_mid))
        ϕZ_mid = zeros(length(∇ϕ1E_mid))
        for i in 1:length(∇ϕ1E_mid)
            E = ∇ϕ1E_mid[i]
            ϕX_mid[i] = E[1]
            ϕZ_mid[i] = E[2]
        end

        ϕX = zeros(length(∇ϕ1E))
        ϕZ = zeros(length(∇ϕ1E))
        for i in 1:length(∇ϕ1E)
            E = ∇ϕ1E[i]
            ϕX[i] = E[1]
            ϕZ[i] = E[2]
        end

        C1_Z = αTosOfFloatLine(S[k][2])


        R = 1/2 * ( ϕX_mid.^2 .+  ϕZ_mid.^2 ) - 9.81.*C1_Z


        #Create ODE step vector
        K = [ ϕX, ϕZ, R ]

        #Step
        S[k+1] = S[k] + nt.*K



    end


    X = [ s[1] for s in S ]
    Y = [ s[2] for s in S ]

    (C2,X,Y)
end
