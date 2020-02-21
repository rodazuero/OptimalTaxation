
cd("C:\\Users\\marya\\Documents\\GitHub\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs")
cd("C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs")

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
#include("NewMyRungeKutta.jl")
include("NewMyRungeKuttaReverse.jl")

#Entrepreneurs Problem
include("NewControlsAlgorithmE.jl")
#include("NewControlsAlgorithmEGrid.jl")
include("States_NewRungeKuttaE.jl")
#include("NewMyRungeKuttaE.jl")
include("NewMyRungeKuttaEReverse.jl")

#The file for the function that computes everything:
include("ProblemFunction.jl")

#Marginal taxes plot
include("marginal_taxes.jl")

#Define values for the model parameters
pa = init_parameters();

#Entrepreneurs Problem
    #Initial boundary conditions (states from the global problem)

    #Define proportion of agents in global problem
    gp     =   0.983

    ue0    =   100.0
    ue0    =   630.0
    μe0    =   0.0 - 1.0e-10
    ye0    =   0.0
    λe0    =   1.0
    le0    =   0.0
    ωe0    =   0.9
    ωe0    =   1.34232

    Nspan = 500
    y_end= [ue0, μe0, ye0, λe0, le0, ωe0, 0.0, 0.0];
    elb= pa.θ_e_b-((1-gp)/pa.he(pa.θ_w_ub,pa.θ_e_ub));
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = eub:-estep:elb;
    solutione = Array{Float64}(undef,Nspan,10);
    my_runge_kuttae_reverse!(solutione,y_end,espan,estep,pa,pa.θ_w_ub)

    #solutione[end,:]
    #using DelimitedFiles
    #writedlm("SolutionNewE.csv",solution,';')

    #For log-normal distribution
    #=xlb= log(pa.θ_w_lb);
    xub = log(pa.θ_w_ub);
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xlb:xstep:xub;
    θspan = exp.(xspan);=#

    θespan = Array{Float64,1}
    θespan = collect(elb:estep:eub)

    controlse = Array{Float64}(undef,Nspan,2);
    recover_controlse!(controlse, pa.θ_w_ub ,θespan, solutione);
    controlse

    #E= DataFrame(controlse)
    #names!(E,[:z, :n] )
    #CSV.write("ControlsEntrepreneurs.csv", E)
    #RungeKuttaE=hcat(DataFrame(solutione),E, DataFrame(thetae=θespan))
    #CSV.write("RungeKuttaNewE.csv", RungeKuttaE)

    #Marginal taxes:
    τ_prime_e = Array{Float64}(undef,Nspan,2);

    marginal_taxese(controlse, θespan, solutione, pa)

    #taxese = DataFrame(τ_prime_e)
    #names!(taxese,[:tau_c, :tau_n] )
    #taxes2e=hcat(DataFrame(theta=θespan),taxese)
    #CSV.write("marginal_taxes_e.csv",taxes2e)

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
    my_runge_kutta_reverse!(solution,y_end,xspan,xstep,pa)

    #using DelimitedFiles
    #writedlm("SolutionNew.csv",solution,';')

    θspan = Array{Float64,1}
    θspan = collect(xlb:xstep:xub)

    controls = Array{Float64}(undef,Nspan,4);
    recover_controls!(controls, θspan, solution);
    controls

    bound_e   = -(pa.indicator*solutione[1,1]^pa.ϕ*he_ub+solutione[1,4]*(e_end*controlse[1,2]^pa.α -(pa.β/(1.0+pa.σ))*controlse[1,1]^(1.0+pa.σ)
                - solutione[1,1])*he_ub - solutione[1,6]*controlse[1,2]*he_ub + solutione[1,2]*controlse[1,2]^pa.α*(1.0-pa.β*controlse[1,1]^pa.σ))

    #C= DataFrame(controls)
    #names!(C,[:z, :n, :l, :p] )
    #RungeKutta=hcat(DataFrame(solution),C, DataFrame(theta=θspan))
    #CSV.write("RungeKuttaNew.csv", RungeKutta)

        #Marginal taxes:
        τ_prime = Array{Float64}(undef,Nspan,3);
        marginal_taxes(controls, θspan, solution, pa)

        #taxes = DataFrame(τ_prime)
        #names!(taxes,[:tau_c, :tau_n, :tau_l] )
        #taxes2=hcat(DataFrame(theta=θspan),taxes)
        #CSV.write("marginal_taxes.csv",taxes2)

        #Utilities:
        utilities_prime = Array{Float64}(undef,Nspan,2);

        for i=1:Nspan
            utilities_prime[i,1] = pa.χ*controls[i,3]^(1.0+pa.ψ)/θspan[i]; #u_w
            utilities_prime[i,2] = controlse[i,2]^pa.α*(1.0-pa.β*controlse[i,1]^pa.σ); #u_e
        end

        #A:
        A_matrix = Array{Float64}(undef,Nspan,3);

        for i=1:Nspan
            A_matrix[i,1] = (pa.indicator*solution[i,1]^pa.ϕ-solution[i,6]*solution[i,1]) + solution[i,4]/pa.he(θspan[i],solution[i,3]); #A
            n_full_info   = ((solution[i,6]*pa.α*solution[i,3])/solution[i,8])^(1.0/(1.0-pa.α));
            A_matrix[i,2] = -(1.0-pa.α)/pa.α*solution[i,8]*n_full_info; #A if z is 0 and n_full_info
            A_matrix[i,3] = controls[i,2]^pa.α*(1.0-τ_prime[i,1]); #u_e in global
        end

        #Plots:
        #graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub,
        #       bound_e,τ_prime,τ_prime_e,"C:\\Users\\marya\\Documents\\GitHub\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs\\Graphs",
        #       utilities_prime,A_matrix)
        graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub,
                bound_e,τ_prime,τ_prime_e,"C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs\\Graphs",
                utilities_prime,A_matrix)
