# Notes and reminders
# 1. g(θ, e) must have several methods allowing both θ and e to be vector or scalars
# 2. When g is the normal distribution, we save a lot of time because the h's are normals
using Distributions
using FastGaussQuadrature
using Cubature #To solve integrals numerically.

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
    uniform_ind::Bool

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
function init_parameters()
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
    ϕ = 0.5;
    G = 0.15;
    #indicator = 0.0; #The Rawlsian case
    uniform_ind = true; #true if we have a uniform distribution
    indicator = 1; #The Utilitarian case

    ## Distributions:
    constant_w_lw = 5.0e-2;
    constant_e_lw = 5.0e-2;
    constant_w_ub = 1.0-0.05;
    constant_e_ub = 1.0-0.05;

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
            θ_e_lb = μ_e-((12.0^0.5)/2.0)*(σ2_e^0.5);
            θ_e_ub = μ_e+((12.0^0.5)/2.0)*(σ2_e^0.5);
            θ_w_lb = μ_w-((12.0^0.5)/2.0)*(σ2_w^0.5);
            θ_w_ub = μ_w+((12.0^0.5)/2.0)*(σ2_w^0.5);

            dist_marginal_w = Uniform(θ_w_lb,θ_w_ub);
            dist_marginal_e = Uniform(θ_e_lb,θ_e_ub);

            cons_mod_dis = 1.0/(1.0-constant_w_lw*constant_e_lw);

    else
        #Log-Normal distribution:
        println("Using Log-Normal Distribution.")

        #Parameters for distributions:
        μ_w = 1.7626;
        μ_e = 1.2528;
        σ2_w = 1.0921; σ_w=σ2_w^0.5;
        σ2_e = 1.1675; σ_e=σ2_e^0.5;
        σ_we = 0.2;
        #σ_we = 0.0;

        #Normal distribution for ln(θ) and log-normal for θw
            dist_marginal_lnθw = Normal(μ_w,σ_w); #This is the distribution for ln(θ_w)
            dist_marginal_lnθe = Normal(μ_e,σ_e); #This is the distribution for ln(θ_e)

            #Lower bounds:
                #θw
                log_θ_w_lb  = quantile(dist_marginal_lnθw,constant_w_lw); #ln(θw)_lb
                θ_w_lb      = exp(log_θ_w_lb);
                #θe
                log_θ_e_lb  = quantile(dist_marginal_lnθe,constant_e_lw); #ln(θe)_lb
                θ_e_lb      = exp(log_θ_e_lb);

            #Upper bounds:
                #θw
                log_θ_w_ub = quantile(dist_marginal_lnθw,constant_w_ub); #ln(θw)_ub
                θ_w_ub     = exp(log_θ_w_ub);
                #θe
                log_θ_e_ub = quantile(dist_marginal_lnθe,constant_e_ub); #ln(θe)_ub
                θ_e_ub     = exp(log_θ_e_ub);

            #Defining the conditional distributions:
                mean_w_given_e(x_e) = μ_w + σ_we*σ_w/σ_e*(x_e-μ_e);
                var_w_given_e = (1.0-σ_we^2.0)*σ2_w;
                dist_w_given_e(x_e) = Normal(mean_w_given_e(x_e),var_w_given_e^0.5);

                mean_e_given_w(x_θ) = μ_e + σ_we*σ_e/σ_w*(x_θ-μ_w);
                var_e_given_w = (1.0-σ_we^2.0)*σ2_e;
                dist_e_given_w(x_θ) = Normal(mean_e_given_w(x_θ),var_e_given_w^0.5);

                cov = σ_we*σ_w*σ_e;
                d   = MvNormal([μ_w, μ_e], [σ2_w cov; cov σ2_e]);
                gg_original(θ,e) = pdf(d,[log(θ),log(e)])/(θ*e);

            #The constant to modify the distributions:
                function integral_1(θ) #We find the value of the first integral, then the double.
                    (val1,err1) = hcubature(x -> gg_original(θ,x[1]), [θ_e_lb],[θ_e_ub]);
                    return val1
                end
                (val,err)    = hcubature(x -> integral_1(x[1]), [θ_w_lb],[θ_w_ub]);
                cons_mod_dis = 1.0/val;
                println(" val_int = ", val, " cons_mod_dis = ", cons_mod_dis)

    end

    #Now we define the functions depending on the distribution we have (they cannot be defined in the conditional above):
    function gg(θ,e)
        if uniform_ind == true
           f = 1.0/((θ_w_ub-θ_w_lb)*(θ_e_ub-θ_e_lb));
        else
           f = gg_original(θ,e)*cons_mod_dis;
        end
    end

    function hw(θ,e)
        if uniform_ind == true
            f = (((e-θ_e_lb)/(θ_e_ub-θ_e_lb))*(1/(θ_w_ub-θ_w_lb))); # h_w(θ)= F_e|w(e|θ) fw(θ)
        elseif σ_we==0.0
            f = 1.0/θ*pdf(dist_marginal_lnθw, log(θ))*(cdf(dist_marginal_lnθe, log(e)) - cdf(dist_marginal_lnθe, log(θ_e_lb)))*cons_mod_dis; # h_w(θ)= f_θw(θ) F_θe(e)
        else
            (val,err) = hcubature(x -> gg(θ,x[1]), [θ_e_lb],[e]);
            f = val;
        end
    end

    function he(θ,e)
        if uniform_ind == true
           f = (((θ-θ_w_lb)/(θ_w_ub-θ_w_lb))*(1/(θ_e_ub-θ_e_lb))); # h_e(e)= F_w|e(θ|e) fe(e)
        elseif σ_we == 0.0
           f = 1.0/e*pdf(dist_marginal_lnθe, log(e))*(cdf(dist_marginal_lnθw, log(θ)) - cdf(dist_marginal_lnθw, log(θ_w_lb)))*cons_mod_dis; # h_e(e)= f_θe(e) F_θw(θ)
        else
           (val,err) = hcubature(x -> gg(x[1],e), [θ_w_lb],[θ]);
           f = val;
        end
    end

    function partialhw_partiale(θ,e)
        f = gg(θ,e); #gg
    end

    function partialhe_partiale(θ,e)
        if uniform_ind == true
            f = 0.0;
        else
            f = ( - 1.0/e*(1.0 + 1.0/(1.0-σ_we^2)*(log(e) - μ_e)/σ2_e )*he(θ,e) +
                    σ_we/(e*σ_e*σ_w*(1.0-σ_we^2))*(pdf(dist_marginal_lnθe, log(e))*(mean_w_given_e(log(e)) -
                    var_w_given_e*pdf(dist_w_given_e(log(e)), ((log(θ) - mean_w_given_e(log(e)))/var_w_given_e))/cdf(dist_w_given_e(log(e)), ((log(θ) - mean_w_given_e(log(e)))/var_w_given_e) ) ) - μ_w*he(θ,e) ) );
        end
    end

    hh(θ,e,p) = hw(θ,e) + p*he(θ,e);

    weights, nodes = gausslegendre(25);

    Param(χ,ψ, κ, ρ, α, δ, γ, β, σ, ς, ϕ, G, indicator, uniform_ind, μ_w, μ_e, σ2_w, σ2_e, σ_we, he, hw, hh, gg, weights, nodes, partialhw_partiale, partialhe_partiale, θ_w_lb, θ_w_ub, θ_e_lb, θ_e_ub, constant_w_lw, constant_e_lw, constant_w_ub, constant_e_ub);
end
