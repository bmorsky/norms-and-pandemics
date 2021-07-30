# SIR time series
using DifferentialEquations, Plots, RCall, SpecialFunctions

tspan = (0.0,50)

# parameters
α = 15
β = 2
γ = 1
μ₁ = 0.7
μ₂ = 0.2
m₁ = 0.7
m₂ = 0.7
σ² = 0.04
ρ = 4
# q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
# q = 0.5*(1 + erf((q̂ - μ)/sqrt(2*σ²)))

δᶜ = 1/6#0.17
k=4

# ODE
function ODE4comp!(du,u,p,t)
    du[1] = -β*u[1]*(u[3]/(1 + (u[5]/δᶜ)^k)+u[4]/(1 + (u[5]/δᶜ)^k)) #/(1 + (u[4]/δᶜ)^k) # Ṡ
    du[2] = -β*u[2]*(u[3]/(1 + (u[5]/δᶜ)^k)+u[4]/(1 + (u[5]/δᶜ)^k)) #/(1 + (u[4]/δᶜ)^k) # Ṡ
    du[3] = - γ*u[3] + β*u[1]*(u[3]/(1 + (u[5]/δᶜ)^k)+u[4]/(1 + (u[5]/δᶜ)^k)) #/(1 + (u[4]/δᶜ)^k) # İ
    du[4] = - γ*u[4] + β*u[2]*(u[3]/(1 + (u[5]/δᶜ)^k)+u[4]/(1 + (u[5]/δᶜ)^k)) #/(1 + (u[4]/δᶜ)^k) # İ
    du[5] = ρ*(0.5*(1 + erf((u[5] - (1 - α*(u[3]+u[4]))*μ₁ - α*(u[3]+u[4])*μ₂)/sqrt(2*σ²))) - u[5]) # ṗ
    du[6] = ρ*(0.5*(1 + erf((u[6] - (1 - α*(u[3]+u[4]))*m₁ - α*(u[3]+u[4])*m₂)/sqrt(2*σ²))) - u[6]) # ṗ
    du[7] = γ*(u[3]+u[4]) # Ṙ
end

u₀ = [0.49;0.49;0.01;0.01;0.0;0.0;0]
prob = ODEProblem(ODE4comp!,u₀,tspan)
sol = solve(prob)
fig = plot(sol,thickness_scaling=2,label = ["S1" "S2" "I1" "I2" "p1" "p2" "R"])
savefig("~/Desktop/2typeshigh.pdf")
