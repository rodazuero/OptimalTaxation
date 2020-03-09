cd("C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\InformalityGlobalAndEntrepreneurs")

using Roots
using NLopt
using Statistics
#using PyPlot
using DataFrames
using CSV
using NLsolve
using Plots
#using DifferentialEquations

#Global Problem
include("Definitions.jl")

pa = init_parameters();

#Entrepreneurs Problem
    #Initial boundary conditions (states from the global problem)

    #Define proportion of agents in global problem
    gp     =   0.993

    ue    =   640.0
    μe    =   0.0 - 1.0e-10
    #μe    =   -0.9
    ye_agg    =   0.0
    λe    =   1.0
    le_agg    =   0.0
    ωfe   =   1.343
    lie_agg   =   0.0
    ωie   =   1.343
    wie   =   1.343
    ϕie   =   0.0


#Construct state object
sse = StateE(ue, μe, ye_agg, λe, le_agg, ωfe, lie_agg, ωie, wie, ϕie);

θe = pa.θ_e_ub
θ = pa.θ_w_ub

    #Recover ditributions
    h_e= pa.he(θ, θe);

    #Defining bounds we use in various cases (limits of z):
    n_full_info = ((sse.λe*pa.α*θe)/sse.ωfe)^(1.0/(1.0-pa.α));
    z_max   = (1.0/pa.β)^(1.0/pa.σ); #Max possible evasion.
    z_lwbar = eps(); #The smallest possible z.
    z_1     = z_max;
    z_2     = -pa.σ*sse.μe*n_full_info^pa.α/(sse.λe*h_e);
    #Keep the lower number:
    z_upbar = min(z_1,z_2);

    #Defining the functions we are using:
    fun_ni(n)   = ((pa.α*θe*n^(pa.α-1.0)-sse.wie)/pa.δ)^(1.0/pa.γ);
    n_opt(z)    = (-sse.λe*z*h_e/(pa.σ*sse.μe))^(1.0/pa.α);
    fun_nz(z,n) = (sse.λe*pa.α*θe*n^(pa.α-1.0)*h_e - sse.ωfe*h_e + pa.α*sse.μe*n^(pa.α-1.0)*(1.0-pa.β*z^pa.σ)+
                   fun_ni(n)^(1.0-pa.γ)/n^(2.0-pa.α)*pa.α*(1.0-pa.α)*θe/(pa.δ*pa.γ)*(sse.λe*pa.δ*fun_ni(n)^pa.γ+sse.ωie-sse.ωfe)*h_e);
    fun_z(z)    = fun_nz(z,n_opt(z));

    #Hamiltonian:
    objective(ze,ne,nie) = pa.indicator*sse.ue^pa.ϕ*h_e + sse.λe*( θe*ne^pa.α - pa.β*ze^(1.0+pa.σ)/(1.0+pa.σ) - pa.δ/(1.0+pa.γ)*nie^(1.0+pa.γ) - sse.ue)*h_e -
                           sse.ωfe*(ne-nei)*h_e + sse.μe*ne^pa.α*(1.0-pa.β*ze^pa.σ) - sse.ωie*nie*h_e;

    #1. The interior solution is:
    #Using bisection method:
    if fun_z(z_lwbar)*fun_z(z_upbar)<0
        potential_z = find_zero(fun_z, (z_lwbar,z_upbar), Bisection());
    else
        error("Bisection: signs equal --> Cannot solve.")
    end
