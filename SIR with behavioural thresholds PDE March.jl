# SIR time series
using DifferentialEquations, Distributions, Plots, RCall, SpecialFunctions

tspan = (0.0,50)

# parameters
α = 1.415
β = 2
γ = 1
ℓ = 4
# q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
# q = 0.5*(1 + erf((q̂ - μ)/sqrt(2*σ²)))

δᶜ = 1/16#0.17
k = 4
M = 20
μ = collect(0.01:0.01:1)
θ = range(0,0.75,length=M)

function calc_gq(u)
    g = 0
    q = 0
    for m = 1:1:M
        δ = 1/(1+exp(-100*(u[2*M+1] - μ[m]/(α+calc_I(u)))))
        g += u[M+m]/(1 + (δ/δᶜ)^k) # I[m]*g[m]
        q += δ*(u[m]+u[M+m])
    end
    return [g, q]
end

function calc_I(u)
    I = 0
    for m = 1:1:M
        I += u[M+m]
    end
    return I
end

function calc_I_beta(u)
    Ibeta = 0
    for m = 1:1:M
        Ibeta += u[M+m]/(1+u[M+2*m])
    end
    return Ibeta
end

# ODE
function PDE(du,u,p,t)
    # u[m] -> S, u[M+m] -> I, u[2*M+m] -> p
    for m = 1:1:M
        du[m] = -β*u[m]*calc_I_beta(u)/(1+u[M+2*m]) # Ṡ[m]
        du[M+m] = -γ*u[M+m] + β*u[m]*calc_I_beta(u)/(1+u[M+2*m]) # İ[m]
        du[2*M+m] = ℓ*(1/(1+exp(-10*(calc_I(u)-0.7+θ[m]))) - u[2*M+m]) # Ḃ[m]
    end
end

S₀ = pdf.(Beta(5,5),range(0,stop=1,length=M))
I₀ = 0.001*S₀
u₀ = vcat(S₀/sum(S₀+I₀),I₀/sum(S₀+I₀),zeros(M))
prob = ODEProblem(PDE,u₀,tspan)
sol = solve(prob,saveat=0.1)
S = sum(sol[1:M,:],dims=1)
I = sum(sol[M+1:2*M,:],dims=1)
R = 1 .- sum(sol[1:2*M,:],dims=1)
p = sum(sol[2*M+1:end,:]/M,dims=1)

fig = plot(sol.t,[S; I; R; p]',thickness_scaling=2,label = ["S" "I" "R" "p"])
savefig("~/Desktop/test3.pdf")
