
cd("C:\\Users\\marya\\Documents\\GitHub\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs")
cd("C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs\\Check_NumericalError")

using Roots
using NLopt
using Statistics
using PyPlot
using DataFrames
using CSV
using NLsolve
#using Plots
#using DifferentialEquations

#Global Problem
include("Definitions2.jl")
include("MJControlsAlgorithmGlobal3.jl")
include("States_NewRungeKutta.jl")
include("NewMyRungeKuttaReverse.jl")
#Entrepreneurs Problem
include("NewControlsAlgorithmE.jl")
include("States_NewRungeKuttaE.jl")
include("NewMyRungeKuttaEReverse.jl")

#Define values for the model parameters
pa = init_parameters();

#Entrepreneurs Problem
    #Define proportion of agents in global problem
    gp     =   0.993
    ue0    =   640.0
    μe0    =   0.0 - 1.0e-10
    ye0    =   0.0
    λe0    =   1.0
    le0    =   0.0
    ωe0    =   1.34297

    Nspan = 500
    y_end= [ue0, μe0, ye0, λe0, le0, ωe0, 0.0, 0.0];
    elb = pa.θ_e_ub - ((1.0-gp)*(pa.θ_e_ub-pa.θ_e_a)*(1.0-pa.constant_w_lw*pa.constant_e_lw));
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = eub:-estep:elb;
    solutione = Array{Float64}(undef,Nspan,10);
    my_runge_kuttae_reverse!(solutione,y_end,espan,estep,pa,pa.θ_w_ub)

#Global Problem (Reverse)
    μ_end = solutione[1,2];
    e_end = elb;
    μ_end = 0.0;
    e_end = pa.θ_e_ub;
    par   = 0.0;
    upper_cons = 0.0
    cons  = ((elb-(pa.θ_e_lb+upper_cons))/((pa.θ_w_ub-pa.θ_w_lb)^(1.0+par)));

    Nspan = 500
    y_end = [μ_end, e_end];
    xlb   = pa.θ_w_lb;
    xub   = pa.θ_w_ub;
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xub:-xstep:xlb;
    solution = Array{Float64}(undef,Nspan,2);
    my_runge_kutta_reverse!(solution,y_end,xspan,xstep,pa, par, cons)

    θspan = Array{Float64,1}
    θspan = collect(xlb:xstep:xub)

    controls = Array{Float64}(undef,Nspan,1);
    recover_controls!(controls, θspan, solution, par, cons);
    controls

    real = Array{Float64}(undef,Nspan,1);

    for i=1:Nspan
        real[i,1] = cons*(θspan[i]-pa.θ_w_lb)^(1.0+par)+(pa.θ_e_lb+upper_cons)
    end

    fig, figura=plt.subplots(1,4)
    fig.suptitle("Figura")
        #μ:
    figura[1].plot(θspan[1:500], solution[1:500,1])
    figura[1].set(ylabel="μ")
        #e:
    figura[2].plot(θspan[1:500], solution[1:500,2])
    figura[2].set(ylabel="e")
        #e real:
    figura[3].plot(θspan[1:500], real[:,1])
    figura[3].set(ylabel="e real")
        #p:
    figura[4].plot(θspan[1:500], controls[1:500,1])
    figura[4].set(ylabel="p")
