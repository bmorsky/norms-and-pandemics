# SIR time series
using DifferentialEquations, Plots, RCall, SpecialFunctions

tspan = (0.0,50)

# parameters
α = 15
β = 2
γ = 1
μ₁ = 0.7
μ₂ = 0.2
σ² = 0.04
ρ = 4
# q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
# q = 0.5*(1 + erf((q̂ - μ)/sqrt(2*σ²)))

δᶜ = 1/6#0.17
k=4

# ODE
function ODE4comp!(du,u,p,t)
    du[1] = -β*u[1]*u[2]*(1-u[4]) #/(1 + (u[4]/δᶜ)^k) # Ṡ
    du[2] = - γ*u[2] + β*u[1]*u[2]*(1-u[4]) #/(1 + (u[4]/δᶜ)^k) # İ
    du[3] = γ*u[2] # Ṙ
    du[4] = ρ*(0.5*(1 + erf((u[4] - (1 - α*u[2])*μ₁ - α*u[2]*μ₂)/sqrt(2*σ²))) - u[4]) # ṗ
end

u₀ = [0.99;0.01;0.0;0]
prob = ODEProblem(ODE4comp!,u₀,tspan)
sol = solve(prob)
fig = plot(sol,thickness_scaling=2,label = ["S" "I" "R" "p"])
savefig("~/Desktop/fig2.pdf")
