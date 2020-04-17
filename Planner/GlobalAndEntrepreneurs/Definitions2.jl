# Notes and reminders
# 1. g(θ, e) must have several methods allowing both θ and e to be vector or scalars
# 2. When g is the normal distribution, we save a lot of time because the h's are normals
using Distributions
using FastGaussQuadrature

mutable struct Param
    ##Workers parameters
    χ::Float64 #Scale parameter for labor supply
    ψ::Float64 #Elasticity of labor supply
    κ::Float64 #Informal supply scale parameter
    ρ::Float64 #Informal supply elasticity

    ## Entrepreneurs parameters
    α::Float64 #Returns to scale/labor elasticity of output
    δ::Float64 #Informal demand scale parameter
    γ::Float64 #Informal demand elasticity
    β::Float64 #Scale parameter for evasion
    σ::Float64 #Elasticity of evasion
    ς::Float64 #self-employed entrepreneurs

    #Planner parameters
    ϕ::Float64 #Concave utilitarian parameter
    G::Float64 #Expenditure needs
    indicator::Float64

    ## Distributions and parameters
    μ_w::Float64
    μ_e::Float64
    σ2_w::Float64
    σ2_e::Float64
    σ_we::Float64
    he::Function
    hw::Function
    hh::Function
    gg::Function
    weights::Array{Float64,1} #añadir tamaño
    nodes::Array{Float64,1}

    #Span of theta
    θ_w_lb::Float64
    θ_w_ub::Float64
    θ_e_lb::Float64
    θ_e_ub::Float64
    constant_w_lw::Float64
    constant_e_lw::Float64
    constant_w_ub::Float64
    constant_e_ub::Float64

end

mutable struct State
    e::Float64 #Threshold function
    uw::Float64 #Utility
    ϕ_e::Float64
    μ::Float64
    λ::Float64
    ω::Float64
end

mutable struct Control
    n::Float64
    p::Float64
    z::Float64
    l::Float64
end

mutable struct StateE
    ue::Float64 #Utility
    μe::Float64
    λe::Float64
    ωe::Float64
end

mutable struct ControlE
    ze::Float64
    ne::Float64
end


#Constructor for the parameters
function init_parameters()
    ##Workers parameters
    χ= 2.0192;
    ψ= 0.4528;
    κ= 0.1021;
    ρ= 0.0912;

    ## Entrepreneurs parameters
    α= 0.73;
    δ= 0.12873;
    γ= 0.7341;
    β= 0.2135;
    σ= 0.1827;
    ς= 10.0;

    #Planner parameters
    ϕ = 0.1;
    G = 0.15;
    indicator = 0.0; #The Rawlsian case
    #indicator = 1; #The Utilitarian case

    ## Distributions:
    μ_w = 1.7626;
    μ_e = 1.2528;
    σ2_w = 1.0921; σ_w=σ2_w^0.5;
    σ2_e = 1.1675; σ_e=σ2_e^0.5;
    σ_we = 0.0;

    #μ_w = 10.0;
    #μ_e = 10.0;
    #σ2_w = 6.0; σ_w=σ2_w^0.5;
    #σ2_e = 6.0; σ_e=σ2_e^0.5;
    #σ_we = 0.0;
    constant_w_lw = 1.0e-2;
    constant_e_lw = 1.0e-2;
    #constant_w_lw = 0.2;
    #constant_e_lw = 0.2;
    #constant_w_ub = 1.0-1.0e-2;
    #constant_e_ub = 1.0-1.0e-2;
    constant_w_ub = 1.0-0.15;
    constant_e_ub = 1.0-0.15;

#Normal distribution for ln(θ) and log-normal for θw
    dist_marginal_w = Normal(μ_w,σ_w); #This is the distribution for ln(θ_w)
    dist_marginal_e = Normal(μ_e,σ_e); #This is the distribution for ln(θ_e)

    #Lower bounds:
        #θw
        log_θ_w_lb  = quantile(dist_marginal_w,constant_w_lw); #ln(θw)_lb
        θ_w_lb      = exp(log_θ_w_lb)
        #θe
        log_θ_e_lb  = quantile(dist_marginal_e,constant_e_lw); #ln(θe)_lb
        θ_e_lb      = exp(log_θ_e_lb)

    #Upper bounds:
        #θw
        log_θ_w_ub = quantile(dist_marginal_w,constant_w_ub) #ln(θw)_ub
        θ_w_ub     = exp(log_θ_w_ub)
        #θe
        log_θ_e_ub = quantile(dist_marginal_e,constant_e_ub) #ln(θe)_ub
        θ_e_ub     = exp(log_θ_e_ub)

    #When distributions are independent:
    hw(θ,e) = 1.0/θ*( pdf(dist_marginal_w, log(θ))*cdf(dist_marginal_e, log(e)) ); # h_w(θ)= f_θw(θ) F_θe(e)
    he(θ,e) = 1.0/e*( pdf(dist_marginal_e, log(e))*cdf(dist_marginal_w, log(θ)) ); # h_e(e)= f_θe(e) F_θw(θ)
    hh(θ,e,p) = hw(θ,e) + p*he(θ,e);
    cov = σ_we*σ_w*σ_e;
    d   = MvNormal([μ_w, μ_e], [σ2_w cov; cov σ2_e]);
    gg(θ,e) = pdf(d,[log(θ),log(e)])/(θ*e);

    weights, nodes = gausslegendre(25);

    Param(χ,ψ, κ, ρ, α, δ, γ, β, σ, ς, ϕ, G, indicator, μ_w, μ_e, σ2_w, σ2_e, σ_we, he, hw, hh, gg, weights, nodes, θ_e_a, θ_w_a, θ_w_lb, θ_w_ub, θ_e_lb, θ_e_ub, constant_w_lw, constant_e_lw);
end
