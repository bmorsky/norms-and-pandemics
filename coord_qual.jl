using RCall, Roots, SpecialFunctions

# parameters
α = 0.01
β = 0.41
γ = 0.1
δ = 0.2
ϵ = 4
κ = 100
η = 0.8

θ = 0.03
c = 0.025

function routhHurwitz(μ,σ²,ι,ℓ,ω,φ)

    # compute equilibria
    eq₁ = [1,0,0,0]
    eq₂ = [1,0,0,1]
    eq₃ = [1,0,0,0.5*(1+c/θ)]
    eq₄ = [1,0,0,1]


    # elements of the Jacobian
    j₁ = φ .+ ι*I/K
    j₂ = φ .+ ι*S/K
    j₃ = ι*I/K
    j₄ = 2*ω*(q .- p).*y.*I
    j₅ = ω*I.*y.^2
    j₆ = ι*y/K
    j₇ = (p .- q).*(ℓ .- 2*ℓ*y .- 2*ω.*(1 .- y).*y)
    j₈ = y.*(1 .- y).*(ℓ .- ω*y)
    j₉ = (q .- p).*f(p̄)
    j₁₀ = 1 .- y.*f(p̄)

    a₃ = j₁ .+ j₇ .+ j₁₀
    a₂ = j₂.*j₃ .+ j₁.*j₇ .+ j₁.*j₁₀ .+ j₇.*j₁₀ .- j₈.*j₉
    a₁ = j₂.*(j₄.*j₆ .+ j₃.*(j₇ .+ j₁₀)) .+ j₁.*(j₇.*j₁₀ .- j₈.*j₉)
    a₀ = j₂.*(j₆.*(j₅.*j₉ .+ j₄.*j₁₀) .+ j₃.*(j₇.*j₁₀ .- j₈.*j₉))

    Δ₄ = a₁.*(a₂.*a₃ .- a₀) .- a₀.*a₃.^2 .> 0
    Δ₃ = a₃ .> 0
    Δ₂ = a₂ .> 0
    Δ₁ = a₁ .> 0
    Δ₀ = a₀ .> 0

    return [length(p), length(p[Δ₄.*Δ₃.*Δ₂.*Δ₁.*Δ₀ .== 1])]
end

# test to see if the population can crash
function crashTest(μ,σ²,ι,ℓ,ω,φ)
    # functions
    soly(x) = (μ.-q .+ sqrt(2*σ²)*erfinv.(2*x .- 1))./(x.-q)
    solỹ(x) = -μ.+q .- sqrt(2*σ²)*erfinv.(2*x .- 1) # (q-p)*soly

    # compute crash equilibria y* and p* within [0,1]
    p = find_zeros(p -> ω*solỹ(p)*soly(p)^2 - (ℓ+ω)*solỹ(p)*soly(p) - ι*soly(p) + ℓ*solỹ(p), 0, q, no_pts=50)
    p = p[0 .<= p .<= 1]
    y = soly(p)

    p = p[0 .<= y .<= 1]
    y = y[0 .<= y .<= 1]

    ℛ₀ = ι./(ω*(q .- p).*y.^2)

    # check if the conditions are met
    if any(x->x<1,ℛ₀)
        return true
    else
        return false
    end
end

function categorize(a,b)
    if ~b
        return  a[2]+1
    else
        return  10+a[2]+1
    end
end

num = 200
ιVℓ = Array{Float64}(undef, num^2, 3)
ιVω = Array{Float64}(undef, num^2, 3)
ιVφ = Array{Float64}(undef, num^2, 3)

ℓVω = Array{Float64}(undef, num^2, 3)
φVℓ = Array{Float64}(undef, num^2, 3)

φVω = Array{Float64}(undef, num^2, 3)

counter=1
for m = 1:1:num
    for n = 1:1:num
        M₁ = 0.2*m/num
        M₂ = m/num
        N₁ = 0.2*n/num
        N₂ = n/num
        ιVℓ[counter,:] = [M₁ N₂ categorize(routhHurwitz(μ,σ²,M₁,N₂,0.1,0.05), crashTest(μ,σ²,M₁,N₂,0.1,0.05))]
        ιVω[counter,:] = [M₁ N₂ categorize(routhHurwitz(μ,σ²,M₁,1,N₂,0.05), crashTest(μ,σ²,M₁,1,N₂,0.05))]
        ιVφ[counter,:] = [M₁ N₁ categorize(routhHurwitz(μ,σ²,M₁,1,0.1,N₁), crashTest(μ,σ²,M₁,1,0.1,N₁))]
        ℓVω[counter,:] = [M₂ N₂ categorize(routhHurwitz(μ,σ²,0.05,M₂,N₂,0.05), crashTest(μ,σ²,0.05,M₂,N₂,0.05))]
        φVℓ[counter,:] = [M₁ N₂ categorize(routhHurwitz(μ,σ²,0.05,N₂,0.1,M₁), crashTest(μ,σ²,0.05,N₂,0.1,M₁))]
        φVω[counter,:] = [M₁ N₂ categorize(routhHurwitz(μ,σ²,0.05,1,N₂,M₁), crashTest(μ,σ²,0.05,1,N₂,M₁))]
        global counter = counter + 1
    end
end

iotaVell = ιVℓ
iotaVomega = ιVω
iotaVphi = ιVφ
ellVomega = ℓVω
phiVell = φVℓ
phiVomega = φVω

@rput iotaVell iotaVomega iotaVphi ellVomega phiVell phiVomega
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

iotaVell <- as.data.frame(iotaVell)
iotaVomega <- as.data.frame(iotaVomega)
iotaVphi <- as.data.frame(iotaVphi)
ellVomega <- as.data.frame(ellVomega)
phiVell <- as.data.frame(phiVell)
phiVomega <- as.data.frame(phiVomega)

cols = c("0"="#2EC4b6", "1"="#F5F5F5", "2"="#011627", "10"="#E71D36", "11"="#FF9F1C")
iotaphicols = c(0,0.05,0.1,0.15,0.2)
iotaphilims = c(0,0.22)
ellomegacols = c(0,0.25,0.5,0.75,1)
ellomegalims = c(0,1.1)

q1 <- ggplot() +
geom_raster(data=iotaVell,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Learning rate, ", '\u2113'))) +
xlab(expression(paste("Inflow rate, ", iota))) +
coord_fixed(ratio = 0.2)
ggsave(q1,filename="coord_iotaVell.png", width = 3.5, height = 3.5)

q2 <- ggplot() +
geom_raster(data=iotaVomega,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Outflow rate, ", omega))) +
xlab(expression(paste("Inflow rate, ", iota))) +
coord_fixed(ratio = 0.2)
ggsave(q2,filename="coord_iotaVomega.png", width = 3.5, height = 3.5)

q3 <- ggplot() +
geom_raster(data=iotaVphi,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
ylab(expression(paste("Resusceptibility rate, ", phi))) +
xlab(expression(paste("Inflow rate, ", iota))) +
coord_fixed(ratio = 1)
ggsave(q3,filename="coord_iotaVphi.png", width = 3.5, height = 3.5)

q4 <- ggplot() +
geom_raster(data=ellVomega,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Outflow rate, ", omega))) +
xlab(expression(paste("Learning rate, ", '\u2113'))) +
coord_fixed(ratio = 1)
ggsave(q4,filename="coord_ellVomega.png", width = 3.5, height = 3.5)

q5 <- ggplot() +
geom_raster(data=phiVell,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Learning rate, ", '\u2113'))) +
xlab(expression(paste("Resusceptibility rate, ", phi))) +
coord_fixed(ratio = 0.2)
ggsave(q5,filename="coord_phiVell.png", width = 3.5, height = 3.5)

q6 <- ggplot() +
geom_raster(data=phiVomega,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Outflow rate, ", omega))) +
xlab(expression(paste("Resusceptibility rate, ", phi))) +
coord_fixed(ratio = 0.2)
ggsave(q6,filename="coord_phiVomega.png", width = 3.5, height = 3.5)

## Function to extract legend
g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

legend <- ggplot() +
geom_raster(data=ellVomega,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="bottom") +
scale_fill_manual(labels=c("0"="cycle","1"="1 stable","2"="2 stable","10"="crash","11"="stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Learning rate, ", '\u2113'))) +
xlab(expression(paste("Inflow rate, ", omega))) +
coord_fixed(ratio = 1) +
guides(fill=guide_legend(title="Nature of equilibria: "))

legend <- g_legend(legend)
ggsave(legend,filename="coord_legend.png", width = 8.5, height = 3.5)
"""
