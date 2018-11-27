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
    weights::Array{Float64,1}
    nodes::Array{Float64,1}

end

mutable struct State
    e::Float64 #Threshold function
    uw::Float64 #Utility
    ϕ_e::Float64
    μ::Float64
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
    α= 0.8; #Garicano et al.
    δ= 0.12873;
    γ= 0.7341;
    β= 0.2135;
    σ= 0.1827;

    #Planner parameters
    ϕ=0.5;
    G=1.0;

    ## Distributions
    #Here I assume normality
    μ_w = 1.7626;
    μ_e = 1.2528;
    σ2_w = 1.0921; σ_w=σ2_w^0.5;
    σ2_e = 1.1675; σ_e=σ2_e^0.5;
    σ_we = 0.2782;
    dist_marginal_w=Normal(μ_w,σ2_w)
    dist_marginal_e=Normal(μ_e,σ2_e)
    mean_w_given_e(e) = μ_w + σ_we*σ_w/σ_e*(e-μ_e)
    var_w_given_e = (1-σ_we^2)σ2_w;
    dist_w_given_e(e)=Normal(mean_w_given_e(e),var_w_given_e);
    mean_e_given_w(θ) = μ_e + σ_we*σ_e/σ_w*(θ-μ_w);
    var_e_given_w= (1-σ_we^2)σ2_e;
    dist_e_given_w(θ)=Normal(mean_e_given_w(θ),var_e_given_w);
    hw(θ,e)=  ( cdf(dist_e_given_w(θ),e) )*( pdf(dist_marginal_w,θ) ); # h_w(θ)= F_e|w(e|θ) fw(θ)
    he(θ,e)=  ( cdf(dist_w_given_e(e),θ) )*( pdf(dist_marginal_e,e) ); # h_e(e)= F_w|e(θ|e) fe(e)
    hh(θ,e,p)= hw(θ,e) + p*he(θ,e);
    d=MvNormal([μ_w, μ_e], [σ2_w σ_we^2; σ_we^2 σ2_e]);
    gg(θ,e) = pdf(d,[θ,e])
    weights, nodes = gausslegendre(20);

    Param(χ,ψ, κ, ρ, α, δ, γ, β, σ, ϕ, G, μ_w, μ_e, σ2_w, σ2_e, σ_we, he, hw, hh, gg, nodes, weights);
end
