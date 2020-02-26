
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

#Define values for the model parameters
pa = init_parameters();

#Entrepreneurs Problem
    #Initial boundary conditions (states from the global problem)

#Global Problem (Reverse)
    gp    = 0.993
    μ_end = -(1-gp)
    e_end = pa.θ_e_ub - ((1-gp)*(pa.θ_e_ub-pa.θ_e_a)*(1.0-pa.constant_w_lw*pa.constant_e_lw));
    he_ub = pa.he(pa.θ_w_ub,e_end)

    Nspan = 500
    y_end= [μ_end, e_end];
    xlb= pa.θ_w_lb;
    xub = pa.θ_w_ub;
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xub:-xstep:xlb;
    solution = Array{Float64}(undef,Nspan,2);
    my_runge_kutta_reverse!(solution,y_end,xspan,xstep,pa)

    θspan = Array{Float64,1}
    θspan = collect(xlb:xstep:xub)

    controls = Array{Float64}(undef,Nspan,1);
    recover_controls!(controls, θspan, solution);
    controls

    fig, figura=plt.subplots(1,3)
    fig.suptitle("Figura")
        #μ:
    figura[1].plot(θspan[1:500], solution[1:500,1])
    figura[1].set(ylabel="μ")
        #e:
    figura[2].plot(θspan[1:500], solution[1:500,2])
    figura[2].set(ylabel="e")
        #p:
    figura[3].plot(θspan[1:500], controls[1:500,1])
    figura[3].set(ylabel="p")
