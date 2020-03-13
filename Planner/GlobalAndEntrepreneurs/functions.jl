cd("C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs")

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
include("Definitions2.jl")
include("MJControlsAlgorithmGlobal3.jl")
include("States_NewRungeKutta.jl")
include("NewMyRungeKuttaReverse.jl")

#Entrepreneurs Problem
include("NewControlsAlgorithmE.jl")
include("States_NewRungeKuttaE.jl")
#include("NewMyRungeKuttaEReverse.jl")
include("TryEntrepreneurs.jl")

#The file for the function that computes everything:
include("ProblemFunction.jl")

#Marginal taxes:
include("marginal_taxes.jl")

#Propositions:
include("Propositions.jl")
include("Integrals.jl")

#Define values for the model parameters
pa = init_parameters();

#Entrepreneurs Problem
    #Initial boundary conditions (states from the global problem)

    #Define proportion of agents in global problem
    gp     =   0.993

    #ue0    =   100.0
    ue0    =   640.0
    #ue0    =   600.0
    μe0    =   0.0 - 1.0e-10
    ye0    =   0.0
    λe0    =   1.0
    le0    =   0.0
    ωe0    =   1.3429705
    ωe0    =   1.35

    Nspan = 500
    y_end= [ue0, μe0, ye0, λe0, le0, ωe0, 0.0, 0.0];
    elb = pa.θ_e_ub - ((1-gp)*(pa.θ_e_ub-pa.θ_e_a)*(1.0-pa.constant_w_lw*pa.constant_e_lw));
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = eub:-estep:elb;
    solutione = Array{Float64}(undef,Nspan,10);
    fill!(solutione,NaN);
    my_runge_kuttae_reverse!(solutione,y_end,espan,estep,pa,pa.θ_w_ub)

    θespan = Array{Float64,1}
    θespan = collect(elb:estep:eub)

    controlse = Array{Float64}(undef,Nspan,2);
    fill!(controlse,NaN);
    recover_controlse!(controlse, pa.θ_w_ub ,θespan, solutione);
    controlse

    #Marginal taxes and taxes path:
    τ_prime_e = Array{Float64}(undef,Nspan,2);
    fill!(τ_prime_e,NaN);
    marginal_taxese(controlse, θespan, solutione, pa);

#Global Problem (Reverse)
    uw_end    = solutione[1,1] #guess
    μ_end     = solutione[1,2]
    e_end     = elb;
    he_ub     = pa.he(pa.θ_w_ub,e_end)
    ϕ_e_end   = -(pa.indicator*solutione[1,1]^pa.ϕ*he_ub+solutione[1,4]*(e_end*controlse[1,2]^pa.α -(pa.β/(1.0+pa.σ))*controlse[1,1]^(1.0+pa.σ)
                - solutione[1,1])*he_ub - solutione[1,6]*controlse[1,2]*he_ub + solutione[1,2]*controlse[1,2]^pa.α*(1.0-pa.β*controlse[1,1]^pa.σ))
    y_agg_end = solutione[1,3]
    λ_end     = 1.0 #guess
    l_agg_end = solutione[1,5]
    ω_end     = solutione[1,6] #guess=#
    l_new_end = 0.0
    y_new_end = 0.0
    ystar_end = 0.0
    lstar_end = 0.0

    Nspan = 500
    y_end= [uw_end, μ_end, e_end, ϕ_e_end, y_agg_end, λ_end, l_agg_end, ω_end, l_new_end, y_new_end, 0, 0 ];
    xlb= pa.θ_w_lb;
    xub = pa.θ_w_ub;
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xub:-xstep:xlb;
    solution = Array{Float64}(undef,Nspan,12);
    fill!(solution,NaN);
    my_runge_kutta_reverse!(solution,y_end,xspan,xstep,pa)

    θspan = Array{Float64,1}
    θspan = collect(xlb:xstep:xub)

    controls = Array{Float64}(undef,Nspan,4);
    fill!(controls,NaN);
    recover_controls!(controls, θspan, solution);
    controls

#Construct state object
#ss = State(e, uw, ϕ_e, μ, λ, ω);
ss = State(solution[268,3], solution[268,1], solution[268,4], solution[268,2], 1.0, ωe0);
θ = θspan[268]

h_e = pa.he(θ, ss.e);
h_w = pa.hw(θ, ss.e);
n_full_info = ((ss.λ*pa.α*ss.e)/ss.ω)^(1.0/(1.0-pa.α));
n_lwbar = 1.0e-10;
n_upbar = n_full_info;
A_cons = (pa.indicator*ss.uw^pa.ϕ-ss.λ*ss.uw)+ss.ϕ_e/h_e;

fun_zagrzero(n)     = ss.e*n^pa.α - ss.ω/ss.λ*n;
lkappa_nzagrzero(n) = pa.σ*pa.β*fun_zagrzero(n)^(pa.σ-1.0)/(1.0-pa.β*fun_zagrzero(n)^pa.σ)*(A_cons +
                      ss.λ*ss.e*n^pa.α - ss.ω*n - ss.λ*fun_zagrzero(n)/pa.σ*(1.0-pa.β/(1.0+pa.σ)*fun_zagrzero(n)^pa.σ));
funct_nagrzero(n)   = pa.α*(-(1.0-pa.α)/pa.α*ss.ω*n + ss.λ/(1.0+pa.σ)*pa.β*fun_zagrzero(n)^(1.0+pa.σ) - A_cons) +
                      lkappa_nzagrzero(n)*(pa.α*ss.e*n^pa.α - ss.ω/ss.λ*n);

funct_nagrzero(n_lwbar)
funct_nagrzero(n_upbar)
plot(funct_nagrzero,n_lwbar,n_upbar)

1e-09*funct_nagrzero(1e-09)
1e-08*funct_nagrzero(1e-08)
1e-07*funct_nagrzero(1e-07)
-pa.α*A_cons
pa.α/n_upbar*(-(1.0-pa.α)/pa.α*ss.ω*n_upbar + ss.λ/(1.0+pa.σ)*pa.β*fun_zagrzero(n_upbar)^(1.0+pa.σ) - A_cons)

pa.α/n_lwbar*(-(1.0-pa.α)/pa.α*ss.ω*n_lwbar + ss.λ/(1.0+pa.σ)*pa.β*fun_zagrzero(n_lwbar)^(1.0+pa.σ) - A_cons)
