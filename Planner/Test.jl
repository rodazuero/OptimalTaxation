cd("C:\\Users\\d-wills\\Dropbox (Uniandes)\\1_5_TaxesAndInformality\\PlannerPblmCodes")

using Roots
using Plots
using DifferentialEquations


include("Definitions.jl")
include("Controls.jl")

#Parameters
pa = init_parameters();
#Prices
ω=0.2;
λ=0.05;
#Current state
θ = exp(pa.μ_w);
e =  exp(pa.μ_e);
uw = 81.0;
μ = 0.1;
ϕ_e= -0.05;
ss=State(e, uw, ϕ_e, μ);
#initialize control objective
ctrl=Array{Control}(undef,1);



 find_controls!(ctrl, θ, λ, ω, ss, pa)
 integrate_dg_de(θ,e,pa)
