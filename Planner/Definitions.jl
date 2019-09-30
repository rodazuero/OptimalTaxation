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
    ϵ::Float64 #self-employed labor supply

    #Planner parameters
    ϕ::Float64 #Concave utilitarian parameter
    G::Float64 #Expenditure needs

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
    ϵ= 0.0;

    #Planner parameters
    ϕ=0.5;
    G=0.15;

    ## Distributions
    #Here I assume normality
    μ_w = 1.7626;
    μ_e = 1.2528;
    #σ2_w = 1.0921; σ_w=σ2_w^0.5;
    σ2_w = 5; σ_w=σ2_w^0.5;
    #σ2_e = 1.1675; σ_e=σ2_e^0.5;
    σ2_e = 5; σ_e=σ2_e^0.5;
    σ_we = 0.0;

    dist_marginal_w=Normal(μ_w,σ_w);
    dist_marginal_e=Normal(μ_e,σ_e);

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

    mean_w_given_e(x_e) = μ_w + σ_we*σ_w/σ_e*(x_e-μ_e);
    var_w_given_e = (1.0-σ_we^2.0)*σ2_w;
    dist_w_given_e(x_e)=Normal(mean_w_given_e(x_e),var_w_given_e^0.5);
    mean_e_given_w(x_θ) = μ_e + σ_we*σ_e/σ_w*(x_θ-μ_w);
    var_e_given_w= (1.0-σ_we^2.0)*σ2_e;
    dist_e_given_w(x_θ)=Normal(mean_e_given_w(x_θ),var_e_given_w^0.5);
    hw(θ,e)= (1/θ) * ( cdf(dist_e_given_w(log(θ)) ,log(e) )*( pdf(dist_marginal_w, log(θ)) ) ); # h_w(θ)= F_e|w(e|θ) fw(θ)
    he(θ,e)= (1/e) * ( cdf(dist_w_given_e(log(e)), log(θ) )*( pdf(dist_marginal_e, log(e)) ) ); # h_e(e)= F_w|e(θ|e) fe(e)
    hh(θ,e,p)= hw(θ,e) + p*he(θ,e);
    cov= σ_we*σ_w*σ_e;
    d=MvNormal([μ_w, μ_e], [σ2_w cov; cov σ2_e]);
    gg(θ,e) = pdf(d,[log(θ),log(e)]) /(θ*e);
    weights, nodes = gausslegendre(25);

    Param(χ,ψ, κ, ρ, α, δ, γ, β, σ, ϵ, ϕ, G, μ_w, μ_e, σ2_w, σ2_e, σ_we, he, hw, hh, gg, nodes, weights, θ_w_lb, θ_w_ub, θ_e_lb, θ_e_ub );
end

#=function distr_hw(θ::Real, e::Real, pa::Param)
    (N,)=size(pa.nodes);
    b = e;
    integral=0.0;
    for i=1:N
        x = (pa.nodes[i] +1.0)*b/2.0;
        integral += pa.gg(θ,x)*pa.weights[i];
    end
    integral = integral*b/2.0;
end


function distr_he(θ::Real, e::Real, pa::Param)
    (N,)=size(pa.nodes);
    b = θ;
    integral=0.0;
    for i=1:N
        x = (pa.nodes[i] +1.0)*b/2.0;
        integral += pa.gg(x,e)*pa.weights[i];
    end
    integral = integral*b/2.0;
end=#
