# SIR time series
using DifferentialEquations, Plots, RCall, SpecialFunctions

tspan = (0.0,20)

# parameters
β = 1/2
μ = 1/2
γ = 1/6
γH = 1/18
fD = 0.5
N = 1e7
δc = 100/N
Dc = 5000/N
k = 2
ϵ = 1/7

# ODE
function ODE4comp!(du,u,p,t)
    du[1] = -β*u[1]*u[3]/(1+(du[5]/δc)^k) # Ṡ
    du[2] = - μ*u[2] + β*u[1]*u[3]/(1+(du[5]/δc)^k) # Ė
    du[3] = μ*u[2] - γ*u[3] # İ
    du[4] = (1-fD)*γ*u[3] # Ṙ
    du[5] = fD*γ*u[3] # Ḋ
end

u₀ = [1e7-1;0;1;0;0]
prob = ODEProblem(ODE4comp!,u₀,tspan)
sol = solve(prob)
fig = plot(sol,thickness_scaling=2,label = ["S" "E" "I" "R" "D"])
savefig("~/Desktop/2typeshigh.pdf")
