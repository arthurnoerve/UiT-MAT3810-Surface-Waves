include("twovector.jl")
include("functions.jl")
include("matrix_elements.jl")
include("boundary.jl")


using ODE
using PyPlot

#Parameters
L = 5; ns = 100; ds = 2L/ns
T = 5; nt = 100; dt = T/nt



#Bottom
#It is described as a parametrical function of the spatial variable s and the time t.

#Returns [X(s,t), Z(x,t)]
function h(s,t)
    X = s
    Z = -4 + sin(π/5*t)exp(-s.^2) - 0.2s
    return [X,Z]
end


bottom,X,Y = boundaryIntegral(L,ns,T,nt,H)


see(100)



function see(n)
    plot(X[n],Y[n]) #Surface
    plot(bottom.X[n,:],bottom.Z[n,:]) #Bottom
end


















# TEST -------------------------------------------------------------------------------



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


k = 1

#Extract
ϕ1E = S[k][3]
∇nϕ2E = ∇nϕ2[k,:]
mean(∇nϕ2E)
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
