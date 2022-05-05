# SEIRB time series
using DifferentialEquations, Distributions, Plots, RCall, SpecialFunctions

# parameters/initial conditions
β = 0.41
δ = 0.2
γ = 0.1
ϵ = 100
pₘ = 0.8
N = 20
f = pdf.(Beta(2,2),range(0,stop=1,length=N))/sum(pdf.(Beta(5,5),range(0,stop=1,length=N)))
Δt = 0.1
T = 1000

function λ(u)
    xpI = 0
    for x = 1:1:N
        xpI += f[x]*u[2*N+x]*(1-pₘ*u[3*N+x])
    end
    return β*xpI
end

function avgI(u)
    y = 0
    for x = 1:1:N
        y += f[x]*u[2*N+x]
    end
    return y
end

function avgp(u)
    y = 0
    for x = 1:1:N
        y += f[x]*u[3*N+x]
    end
    return y
end

function BR(u,x)
    xp = 0
    for x = 1:1:N
        xp += f[x]*(1-pₘ*u[3*N+x])
    end
    return 1/(1+exp(-100*(avgp(u)+20*avgI(u)-1))) #1/(1+exp(-100*(x/N+xp+50*λ(u)-2.5)))
end

# PDE
function PDE(du,u,p,t)
    # u[x] = S, u[N+x] = E, u[2N+x] = I, u[3N+x] = p
    for x = 1:1:N
        du[x] = -λ(u)*(1-pₘ*u[3*N+x])*u[x] # Ṡ(x)
        du[N+x] = λ(u)*(1-pₘ*u[3*N+x])*u[x] - δ*u[N+x] #Ė(x)
        du[2*N+x] = δ*u[N+x] - γ*u[2*N+x] #İ(x)
        du[3*N+x] = ϵ*(BR(u,x) - u[3*N+x]) #Ḃ(x)
    end
end

u₀ = vcat(0.99*ones(N), zeros(N), 0.01*ones(N), zeros(N))
prob = ODEProblem(PDE,u₀,(0.0,T))
sol = solve(prob,saveat=Δt)

fig = plot(sol.t,sol[3*N:4*N,:]',thickness_scaling=2)
savefig("~/Desktop/inf_eps1.pdf")

S = zeros(Int(T/Δt+1))
E = zeros(Int(T/Δt+1))
I = zeros(Int(T/Δt+1))
B = zeros(Int(T/Δt+1))
for t = 1:1:Int(T/Δt+1)
    S[t] = sum(f.*sol[1:N,t])
    E[t] = sum(f.*sol[N+1:2*N,t])
    I[t] = sum(f.*sol[2*N+1:3*N,t])
    B[t] = sum(f.*sol[3*N+1:4*N,t])
end

fig = plot(sol.t,[S E I B],thickness_scaling=2,label = ["S" "E" "I" "B"])
savefig("~/Desktop/test_eps1beta4.pdf")
#
# out = zeros(N)
# for θ = 1:1:N
#     out .+= f[θ]./(1 .+ exp.(-100*(range(0,stop=1,length=N) .+(θ-1)/N .+1 .-1.5)))
# end
#
# fig = plot(range(0,stop=1,length=N),out,thickness_scaling=2)
# savefig("~/Desktop/test2.pdf")
