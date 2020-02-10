
cd("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Reverse2")
cd("C:\\Users\\mariagon\\Dropbox\\Reverse2")

using Roots
using NLopt
using Statistics
#using PyPlot
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
include("ControlsAlgorithmE.jl")
include("States_NewRungeKuttaE.jl")
#include("NewMyRungeKuttaE.jl")
include("NewMyRungeKuttaEReverse.jl")

include("ProblemFunction.jl")



#Marginal taxes plot
include("marginal_taxes.jl")


#Define values for the model parameters
pa = init_parameters();











#Entrepreneurs Problem
    #Initial boundary conditions (states from the global problem)

    #Define proportion of agents in global problem
    gp     =   0.9

    ue0    =   1000.0
    μe0    =   0.0
    ye0    =   0.0
    λe0    =   1.0
    le0    =   0.0
    ωe0    =   2.75

    Nspan = 500
    y_end= [ue0, μe0, ye0, λe0, le0, ωe0, 0.0, 0.0];
    elb= pa.θ_e_b-((1-gp)/pa.he(pa.θ_w_ub,pa.θ_e_ub));
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = eub:-estep:elb;
    solutione = Array{Float64}(undef,Nspan,10);
    my_runge_kuttae_reverse!(solutione,y_end,espan,estep,pa,pa.θ_w_ub)

    solutione[end,:]
    using DelimitedFiles
    writedlm("SolutionNewE.csv",solution,';')

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

    E= DataFrame(controlse)
    names!(E,[:z, :n] )
    CSV.write("ControlsEntrepreneurs.csv", E)
    RungeKuttaE=hcat(DataFrame(solutione),E, DataFrame(thetae=θespan))
    CSV.write("RungeKuttaNewE.csv", RungeKuttaE)

    #Plots:
    graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub, bound_e,τ_prime,τ_prime_e)

        #States
        estados=subplot(321)
        suptitle("Optimal states")
        estados=plot(θespan[1:500], solutione[1:500,1])
        ylabel("u_e")
        estados=subplot(322)
        estados[1,2]=plot(θespan, solutione[:,2])
        plot(θespan,repeat([0],500),color="lime")
        ylabel("μe")
        estados=subplot(323)
        estados[2,1]=plot(θespan, solutione[:,3])
        plot(θespan,repeat([pa.θ_e_ub],500),color="lime")
        ylabel("Ye")
        estados=subplot(324)
        estados[2,2]=plot(θespan, solutione[:,4])
        #plot(θspan,repeat([bound_e],500),color="lime")
        ylabel("λe")
        estados=subplot(325)
        estados[3,1]=plot(θespan, solutione[:,5])
        plot(θespan,repeat([0],500),color="lime")
        ylabel("Le")
        estados=subplot(326)
        estados[3,2]=plot(θespan, solutione[:,6])
        ylabel("ωe")

        #Auxiliar states
        estados=subplot(121)
        estados[1,2]=plot(θespan, solutione[:,9])
        plot(θespan,repeat([0],500),color="lime")
        ylabel("Le*")
        estados=subplot(122)
        estados[1,2]=plot(θespan, solutione[:,10])
        plot(θespan,repeat([0.15],500),color="lime")
        ylabel("Ye*")

        #Controls
        ctrl=subplot(121)
        suptitle("Optimal controls")
        ctrl=plot(θespan, controlse[:,1])
        ylabel("z")
        ctrl=subplot(122)
        ctrl[1,2]=plot(θespan, controlse[:,2])
        ylabel("n")

        #Marginal taxes:
        τ_prime_e = Array{Float64}(undef,Nspan,2);

        marginal_taxese(controlse, θespan, solutione, pa)

        taxese = DataFrame(τ_prime_e)
        names!(taxese,[:tau_c, :tau_n] )
        taxes2e=hcat(DataFrame(theta=θespan),taxese)
        CSV.write("marginal_taxes_e.csv",taxes2e)

        #Marginal taxes plot:
        plot_margtax = subplot(121)
        suptitle("Marginal taxes")
        plot_margtax = plot(θespan[1:500], τ_prime_e[1:500,1])
        ylabel("τ_c′")
        plot_margtax = subplot(122)
        plot_margtax[1,2] = plot(θespan[1:500], τ_prime_e[1:500,2])
        ylabel("τ_n′")



#Global Problem (Reverse)

    include("NewMyRungeKuttaReverse.jl")


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

    using DelimitedFiles
    writedlm("SolutionNew.csv",solution,';')


    θspan = Array{Float64,1}
    θspan = collect(xlb:xstep:xub)

    controls = Array{Float64}(undef,Nspan,4);
    recover_controls!(controls, θspan, solution);
    controls

    bound_e   = -(pa.indicator*solutione[1,1]^pa.ϕ*he_ub+solutione[1,4]*(e_end*controlse[1,2]^pa.α -(pa.β/(1.0+pa.σ))*controlse[1,1]^(1.0+pa.σ)
                - solutione[1,1])*he_ub - solutione[1,6]*controlse[1,2]*he_ub + solutione[1,2]*controlse[1,2]^pa.α*(1.0-pa.β*controlse[1,1]^pa.σ))


    C= DataFrame(controls)
    names!(C,[:z, :n, :l, :p] )
    RungeKutta=hcat(DataFrame(solution),C, DataFrame(theta=θspan))
    CSV.write("RungeKuttaNew.csv", RungeKutta)

    #Graphs Global problem

        #States
        estados=subplot(421)
        suptitle("Optimal states")
        estados=plot(θspan[1:500], solution[1:500,1])
        ylabel("u_w")
        estados=subplot(422)
        estados[1,2]=plot(θspan, solution[:,2])
        plot(θspan,repeat([0],500),color="lime")
        ylabel("μ")
        estados=subplot(423)
        estados[2,1]=plot(θspan, solution[:,3])
        plot(θspan,repeat([pa.θ_e_ub],500),color="lime")
        ylabel("e")
        estados=subplot(424)
        estados[2,2]=plot(θspan, solution[:,4])
        #plot(θspan,repeat([bound_e],500),color="lime")
        plot(θspan,repeat([bound_e],500),color="lime")
        ylabel("φe")
        estados=subplot(425)
        estados[3,1]=plot(θspan, solution[:,5])
        ylabel("Y")
        estados=subplot(426)
        estados[3,2]=plot(θspan, solution[:,6])
        ylabel("λ")
        estados=subplot(427)
        estados[4,1]=plot(θspan, solution[:,7])
        plot(θspan,repeat([0],500),color="lime")
        ylabel("L")
        estados=subplot(428)
        estados[4,2]=plot(θspan, solution[:,8])
        ylabel("ω")

        #Auxiliar states
        estados=subplot(121)
        estados[1,2]=plot(θspan, solution[:,11])
        plot(θspan,repeat([0],500),color="lime")
        ylabel("L*")
        estados=subplot(122)
        estados[1,2]=plot(θspan, solution[:,12])
        plot(θspan,repeat([0.15],500),color="lime")
        ylabel("Y*")

        #Controls
        ctrl=subplot(221)
        suptitle("Optimal controls")
        ctrl=plot(θspan, controls[:,1])
        ylabel("z")
        ctrl=subplot(222)
        ctrl[1,2]=plot(θspan, controls[:,2])
        ylabel("n")
        ctrl=subplot(223)
        ctrl[2,1]=plot(θspan, controls[:,3])
        ylabel("l")
        ctrl=subplot(224)
        ctrl[2,2]=plot(θspan, controls[:,4])
        ylabel("p")

        #Marginal taxes:
        τ_prime = Array{Float64}(undef,Nspan,3);
        marginal_taxes(controls, θspan, solution, pa)

        taxes = DataFrame(τ_prime)
        names!(taxes,[:tau_c, :tau_n, :tau_l] )
        taxes2=hcat(DataFrame(theta=θspan),taxes)
        CSV.write("marginal_taxes.csv",taxes2)

        #Marginal taxes plot:
        plot_margtax = subplot(131)
        suptitle("Marginal taxes")
        plot_margtax = plot(θspan[1:500], τ_prime[1:500,1])
        ylabel("τ_c′")
        plot_margtax = subplot(132)
        plot_margtax[1,2] = plot(θspan[1:500], τ_prime[1:500,2])
        ylabel("τ_n′")
        plot_margtax = subplot(133)
        plot_margtax[1,3] = plot(θspan[1:500], τ_prime[1:500,3])
        ylabel("τ_l′")
