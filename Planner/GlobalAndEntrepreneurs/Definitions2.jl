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

    #Derivatives for states:
    partialhw_partiale::Function
    partialhe_partiale::Function

    #Span of theta
    θ_e_a::Float64
    θ_w_a::Float64
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
function init_parameters(uniform_ind::Bool)
    #Workers parameters
    χ= 2.0192;
    ψ= 0.4528;
    κ= 0.1021;
    ρ= 0.0912;

    #Entrepreneurs parameters
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
    constant_w_lw = 1.0e-2;
    constant_e_lw = 1.0e-2;
    constant_w_ub = 1.0-0.15;
    constant_e_ub = 1.0-0.15;

    #Getting the distributions depending if we have a log-normal or a uniform:
    #First we define the parameters and bounds and then the functions we are using:
    if uniform_ind == true
        #Uniform distribution:
        println("Using Uniform Distribution.")

            #Parameters for distributions:
            μ_w = 10.0;
            μ_e = 10.0;
            σ2_w = 6.0; σ_w=σ2_w^0.5;
            σ2_e = 6.0; σ_e=σ2_e^0.5;
            σ_we = 0.0;

            #Defining the bounds:
            θ_e_a  = μ_e-((12.0^0.5)/2.0)*(σ2_e^0.5);
            θ_e_ub = μ_e+((12.0^0.5)/2.0)*(σ2_e^0.5);
            θ_w_a  = μ_w-((12.0^0.5)/2.0)*(σ2_w^0.5);
            θ_w_ub = μ_w+((12.0^0.5)/2.0)*(σ2_w^0.5);

            dist_marginal_w = Uniform(θ_w_a,θ_w_ub);
            dist_marginal_e = Uniform(θ_e_a,θ_e_ub);
            θ_w_lb = quantile(dist_marginal_w,constant_w_lw);
            θ_e_lb = quantile(dist_marginal_e,constant_e_lw);

            cons_mod_dis = 1.0/(1.0-constant_w_lw*constant_e_lw);
    else
        #Log-Normal distribution:
        println("Using Log-Normal Distribution.")

        #We do not have bounds in these distribution, so we put them in zero:
        θ_w_a = 0.0;
        θ_e_a = 0.0;

        #Parameters for distributions:
        μ_w = 1.7626;
        μ_e = 1.2528;
        σ2_w = 1.0921; σ_w=σ2_w^0.5;
        σ2_e = 1.1675; σ_e=σ2_e^0.5;
        σ_we = 0.2;

        #Normal distribution for ln(θ) and log-normal for θw
            dist_marginal_lnθw = Normal(μ_w,σ_w); #This is the distribution for ln(θ_w)
            dist_marginal_lnθe = Normal(μ_e,σ_e); #This is the distribution for ln(θ_e)

            #Lower bounds:
                #θw
                log_θ_w_lb  = quantile(dist_marginal_lnθw,constant_w_lw); #ln(θw)_lb
                θ_w_lb      = exp(log_θ_w_lb)
                #θe
                log_θ_e_lb  = quantile(dist_marginal_lnθe,constant_e_lw); #ln(θe)_lb
                θ_e_lb      = exp(log_θ_e_lb)

            #Upper bounds:
                #θw
                log_θ_w_ub = quantile(dist_marginal_lnθw,constant_w_ub) #ln(θw)_ub
                θ_w_ub     = exp(log_θ_w_ub)
                #θe
                log_θ_e_ub = quantile(dist_marginal_lnθe,constant_e_ub) #ln(θe)_ub
                θ_e_ub     = exp(log_θ_e_ub)

            #Defining the conditional distributions:
                mean_w_given_e(x_e) = μ_w + σ_we*σ_w/σ_e*(x_e-μ_e);
                var_w_given_e = (1.0-σ_we^2.0)*σ2_w;
                dist_w_given_e(x_e) = Normal(mean_w_given_e(x_e),var_w_given_e^0.5);

                mean_e_given_w(x_θ) = μ_e + σ_we*σ_e/σ_w*(x_θ-μ_w);
                var_e_given_w = (1.0-σ_we^2.0)*σ2_e;
                dist_e_given_w(x_θ) = Normal(mean_e_given_w(x_θ),var_e_given_w^0.5);

                cov     = σ_we*σ_w*σ_e;
                d       = MvNormal([μ_w, μ_e], [σ2_w cov; cov σ2_e]);
    end

    #Now we define the functions depending on the distribution we have (they cannot be defined in the conditional above):
    function hw(θ,e)
        if uniform_ind == true
            f = (((e-θ_e_a)/(θ_e_ub-θ_e_a))*(1/(θ_w_ub-θ_w_a)))*cons_mod_dis; # h_w(θ)= F_e|w(e|θ) fw(θ)
        elseif σ_we == 0.0
            f = 1.0/θ*( pdf(dist_marginal_lnθw, log(θ))*cdf(dist_marginal_lnθe, log(e)) ); # h_w(θ)= f_θw(θ) F_θe(e)
        else
            f = 1.0/θ*( pdf(dist_marginal_lnθw, log(θ))*cdf(dist_e_given_w(log(θ)), log(e)) ); # h_w(θ)= f_θw(θ) F_θe|θ(e|θ)
        end

        return f
    end

    function he(θ,e)
           if uniform_ind == true
               f = (((θ-θ_w_a)/(θ_w_ub-θ_w_a))*(1/(θ_e_ub-θ_e_a)))*cons_mod_dis; # h_e(e)= F_w|e(θ|e) fe(e)
           elseif σ_we == 0.0
               f = 1.0/e*( pdf(dist_marginal_lnθe, log(e))*cdf(dist_marginal_lnθw, log(θ)) ); # h_e(e)= f_θe(e) F_θw(θ)
           else
               f = 1.0/e*( pdf(dist_marginal_lnθe, log(e))*cdf(dist_w_given_e(log(e)), log(θ)) ); # h_e(e)= f_θe(e) F_θw|e(w|e)
          end

          return f
    end

    function gg(θ,e)
           if uniform_dist == true
               f = (1.0/((θ_w_ub-θ_w_a)*(θ_e_ub-θ_e_a)))*cons_mod_dis;
           else
               f = pdf(d,[log(θ),log(e)])/(θ*e);
          end

          return f
    end

    function partialhw_partiale(θ,e)
           if uniform_dist == true
               f = (1.0/((θ_w_ub-θ_w_a)*(θ_e_ub-θ_e_a)))*cons_mod_dis;
           elseif σ_we == 0.0
               f = pdf(dist_marginal_lnθe, log(e))/e*pdf(dist_marginal_lnθw, log(θ))/θ; #f_θw(θ) f_θe(e)
           else
               f = pdf(dist_marginal_lnθw, log(θ))/θ*pdf(d,[log(θ),log(e)])/(θ*e)/pdf(dist_marginal_lnθe, log(e))/e; #f_θw(θ) f_θe(e)
          end

          return f
    end

    function partialhe_partiale(θ,e)
        if uniform_dist == true
           f = 0.0;
        elseif σ_we == 0.0
           f = -pdf(dist_marginal_lnθe, log(e))/e*(1.0+(log(e)-pa.μ_e)/pa.σ2_e)*cdf(dist_marginal_lnθw, log(θ));
        else
            f = -pdf(dist_marginal_lnθe, log(e))/e*(1.0+(log(e)-pa.μ_e)/pa.σ2_e)*cdf(dist_e_given_w(log(θ)), log(e));
        end

        return f
    end

    hh(θ,e,p) = hw(θ,e) + p*he(θ,e);

    weights, nodes = gausslegendre(25);

    Param(χ,ψ, κ, ρ, α, δ, γ, β, σ, ς, ϕ, G, indicator, μ_w, μ_e, σ2_w, σ2_e, σ_we, he, hw, hh, gg, weights, nodes, partialhw_partiale, partialhe_partiale, θ_e_a, θ_w_a, θ_w_lb, θ_w_ub, θ_e_lb, θ_e_ub, constant_w_lw, constant_e_lw, constant_w_ub, constant_e_ub);
end
