# Notes and reminders
# 1. g(θ, e) must have several methods allowing both θ and e to be vector or scalars
# 2. When g is the normal distribution, we save a lot of time because the h's are normals
using Distributions
using FastGaussQuadrature

struct EconomicParameters
    ## Workers' parameters
    χ::Float64 		# Scale parameter for labor supply
    ψ::Float64 		# Elasticity of labor supply
    κ::Float64 		# Informal supply scale parameter
    ρ::Float64 		# Informal supply elasticity
    ## Entrepreneurs' parameters
    α::Float64 		# Returns to scale/labor elasticity of output
    δ::Float64 		# Informal demand scale parameter
    γ::Float64 		# Informal demand elasticity
    β::Float64	 	# Scale parameter for evasion
    σ::Float64 		# Elasticity of evasion
    ς::Float64 		# Entrepreneurs' self-provision of labor hours
    ## Planner's parameters
    ϕ::Float64 		# Concave utilitarian parameter
    G::Float64 		# Expenditure needs
    utilit::Bool 	# Indicator of utilitarian planner
end

struct DistributionParameters
    # Parameters for distributions
    μ_w::Float64	# Mean of worker's ability
    μ_e::Float64	# Mean of entrepreneur's ability
    σ2_w::Float64	# Variance of worker's ability
    σ2_e::Float64	# Varaince of entrepreneur's ability
    σ_we::Float64	# Covariance of abilities.
    # Span of theta
    θ_w_l::Float64	# Lower bound of worker's ability
    θ_w_u::Float64	# Upper bound of worker's ability
    θ_e_l::Float64	# Lower bound of entrepreneurs's ability
    θ_e_u::Float64	# Upper bound of entrepreneurs's ability
    # Type of distribution
    uniform::Bool	# Indicator of uniform distributions.
end

struct ComputationParameters
	gridsize::Int64
	egridsize::Int64
end

struct ModelDistributions
    he::Function
    hw::Function
    gg::Function
    partial_he_e::Function
    function ModelDistributions(distpar::DistributionParameters)
    	if distpar.uniform # Uniform case
    		hw(θvar,evar)	= (evar-distpar.θ_e_l)/(distpar.θ_e_u-distpar.θ_e_l)/(distpar.θ_w_u-distpar.θ_w_l) # h_w(θ,e) = F_e|w(e|θ)*f_w(θ)
    		he(θvar,evar)	= (θvar-distpar.θ_w_l)/(distpar.θ_w_u-distpar.θ_w_l)/(distpar.θ_e_u-distpar.θ_e_l) # h_e(θ,e) = F_w|e(θ|e)*f_e(e)
    		gg(θvar,evar)	= 1.0/(distpar.θ_w_u-distpar.θ_w_l)/(distpar.θ_e_u-distpar.θ_e_l)	# g(θ,e) = f_w|e(θ|e)*f_e(e) = f_e|w(e|θ)*f_w(θ)
    		partial_he_e(θvar,evar)	= 0.0
    	else # Log-Normal case
    		hw(θvar,evar)	= NaN # h_w(θ,e) = F_e|w(e|θ)*f_w(θ)
    		he(θvar,evar)	= NaN # h_e(θ,e) = F_w|e(θ|e)*f_e(e)
    		gg(θvar,evar)	= NaN	# g(θ,e) = f_w|e(θ|e)*f_e(e) = f_e|w(e|θ)*f_w(θ)
    		partial_he_e(θvar,evar)	= NaN	# ∂h_e/∂e
    	end
    	new(he, hw, gg, partial_he_e)
    end
end

struct OTEmodel
	ecopar::EconomicParameters
	compar::ComputationParameters
	dispar::DistributionParameters
	distri::ModelDistributions
	states::Array{Float64,2}
	controls::Array{Float64,2}
	prices::Array{Float64,1}	# Placeholder for constant states
	function OTEmodel(ecopar::EconomicParameters, compar::ComputationParameters, dispar::DistributionParameters)
		distributions=ModelDistributions(dispar)
		states=Array{Float64}(undef,6, compar.gridsize)
		controls=Array{Float64}(undef,6, compar.gridsize)
		prices=Array{Float64}(undef,2)
		new(ecopar, compar, dispar, ditributions,states,controls,prices)
	end	
end

mutable struct SolverParameters
    # Lower tail cuts
    cut_w::Float64	# Mass of lower tail cut on workers'ability dist.
    cut_e::Float64	# Mass of lower tail cut on entrepreneurs'ability dist.
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
    #μ_w = 1.7626;
    #μ_e = 1.2528;
    #σ2_w = 1.0921; σ_w=σ2_w^0.5;
    #σ2_e = 1.1675; σ_e=σ2_e^0.5;
    #σ_we = 0.0;

    μ_w = 10.0;
    μ_e = 10.0;
    σ2_w = 6.0; σ_w=σ2_w^0.5;
    σ2_e = 6.0; σ_e=σ2_e^0.5;
    σ_we = 0.0;
    constant_w_lw = 1.0e-2;
    constant_e_lw = 1.0e-2;

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

    cons_mod_dis = 1.0/(1.0-constant_w_lw*constant_e_lw);
    hw(θvar,e)   = (((e-θ_e_a)/(θ_e_ub-θ_e_a))*(1/(θ_w_ub-θ_w_a)))*cons_mod_dis; # h_w(θ)= F_e|w(e|θ) fw(θ)
    he(θ,e)   = (((θ-θ_w_a)/(θ_w_ub-θ_w_a))*(1/(θ_e_ub-θ_e_a)))*cons_mod_dis; # h_e(e)= F_w|e(θ|e) fe(e)
    hh(θ,e,p) = hw(θ,e) + p*he(θ,e);

    gg(θ,e) = (1.0/((pa.θ_w_ub-pa.θ_w_a)*(pa.θ_e_ub-pa.θ_e_a)))*cons_mod_dis; #For uniform case.

    weights, nodes = gausslegendre(25);

    Param(χ,ψ, κ, ρ, α, δ, γ, β, σ, ς, ϕ, G, indicator, μ_w, μ_e, σ2_w, σ2_e, σ_we, he, hw, hh, gg, weights, nodes, θ_e_a, θ_w_a, θ_w_lb, θ_w_ub, θ_e_lb, θ_e_ub, constant_w_lw, constant_e_lw);
end
