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
    θ_e_a::Float64
    θ_w_a::Float64
    θ_w_lb::Float64
    θ_w_ub::Float64
    θ_e_lb::Float64
    θ_e_ub::Float64
    constant_w_lw::Float64
    constant_e_lw::Float64

end

mutable struct State
    e::Real #Threshold function
    uw::Real #Utility
    ϕ_e::Real
    μ::Real
    λ::Real
    ω::Real
end

mutable struct Control
    n::Real
    p::Real
    z::Real
    l::Real
end

mutable struct StateE
    ue::Real #Utility
    μe::Real
    λe::Real
    ωe::Real
end

mutable struct ControlE
    ze::Real
    ne::Real
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
    constant_w_ub = 1.0-1.0e-2;
    constant_e_ub = 1.0-1.0e-2;

    #Log normal distribution
        dist_marginal_w = Normal(μ_w,σ_w); #This is the distribution for ln(θ_w)
        dist_marginal_e = Normal(μ_e,σ_e); #This is the distribution for ln(θ_e)

        #Lower bounds:
            #θw
            θ_w_x  = quantile(dist_marginal_w,constant_w_lw);
            θ_w_lb = exp(θ_w_x)
            #θe
            θ_e_x  = quantile(dist_marginal_e,constant_e_lw);
            θ_e_lb = exp(θ_e_x)

        #Upper bounds:
            #θw
            θ_w_x  = quantile(dist_marginal_w,constant_w_ub)
            θ_w_lb = exp(θ_w_x)
            #θe
            θ_e_x  = quantile(dist_marginal_e,constant_e_lw)
            θ_e_lb = exp(θ_e_x)

        #theta_e_lb
        quantile_theta_e_lb(k) = cdf(dist_marginal_e,k) - 0.3
        x_e_lb = find_zero(quantile_theta_e_lb, (-100.0,100.0))
        θ_e_lb= exp(x_e_lb)
        #theta_e_ub
        quantile_theta_e_ub(k) = cdf(dist_marginal_e,k) - 0.7
        x_e_ub = find_zero(quantile_theta_e_ub, (-100.0,100.0))
        θ_e_ub= exp(x_e_ub)

#Uniform distribution
    θ_e_a  = μ_e-((12.0^0.5)/2)*(σ2_e^0.5);
    θ_e_ub = μ_e+((12.0^0.5)/2)*(σ2_e^0.5);

    #theta_w_a
    θ_w_a  = μ_w-((12.0^0.5)/2)*(σ2_w^0.5);
    θ_w_ub = μ_w+((12.0^0.5)/2)*(σ2_w^0.5);

    dist_marginal_w = Uniform(θ_w_a,θ_w_ub);
    dist_marginal_e = Uniform(θ_e_a,θ_e_ub);

    #theta_w_lb
    θ_w_lb = quantile(dist_marginal_w,constant_w_lw);

    #theta_e_lb
    θ_e_lb = quantile(dist_marginal_e,constant_e_lw);



    #Log normal distribution
        dist_marginal_w=Normal(μ_w,σ_w); #This is the distribution for ln(θ_w)
        dist_marginal_e=Normal(μ_e,σ_e); #This is the distribution for ln(θ_e)

        #theta_w_lb
        quantile_theta_w_lb(k) = cdf(dist_marginal_w,k) - 0.3
        x_w_lb = find_zero(quantile_theta_w_lb, (-100.0,100.0))
        θ_w_lb= exp(x_w_lb)
        #theta_w_ub
        quantile_theta_w_ub(k) = cdf(dist_marginal_w,k) - 0.7
        x_w_ub = find_zero(quantile_theta_w_ub, (-100.0,100.0))
        θ_w_ub= exp(x_w_ub)


        #theta_e_lb
        quantile_theta_e_lb(k) = cdf(dist_marginal_e,k) - 0.3
        x_e_lb = find_zero(quantile_theta_e_lb, (-100.0,100.0))
        θ_e_lb= exp(x_e_lb)
        #theta_e_ub
        quantile_theta_e_ub(k) = cdf(dist_marginal_e,k) - 0.7
        x_e_ub = find_zero(quantile_theta_e_ub, (-100.0,100.0))
        θ_e_ub= exp(x_e_ub)



    cons_mod_dis = 1.0/(1.0-constant_w_lw*constant_e_lw);
    hw(θ,e)   = (((e-θ_e_a)/(θ_e_ub-θ_e_a))*(1/(θ_w_ub-θ_w_a)))*cons_mod_dis; # h_w(θ)= F_e|w(e|θ) fw(θ)
    he(θ,e)   = (((θ-θ_w_a)/(θ_w_ub-θ_w_a))*(1/(θ_e_ub-θ_e_a)))*cons_mod_dis; # h_e(e)= F_w|e(θ|e) fe(e)
    hh(θ,e,p) = hw(θ,e) + p*he(θ,e);

    gg(θ,e) = (1.0/((pa.θ_w_ub-pa.θ_w_a)*(pa.θ_e_ub-pa.θ_e_a)))*cons_mod_dis; #For uniform case.

    weights, nodes = gausslegendre(25);

    Param(χ,ψ, κ, ρ, α, δ, γ, β, σ, ς, ϕ, G, indicator, μ_w, μ_e, σ2_w, σ2_e, σ_we, he, hw, hh, gg, weights, nodes, θ_e_a, θ_w_a, θ_w_lb, θ_w_ub, θ_e_lb, θ_e_ub, constant_w_lw, constant_e_lw);
end
