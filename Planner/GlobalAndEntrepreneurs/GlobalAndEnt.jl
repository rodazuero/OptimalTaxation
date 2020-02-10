
cd("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\GlobalAndEntrepreneursU")
cd("C:\\Users\\mariagon\\Dropbox\\GlobalAndEntrepreneursU")

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
include("ControlsAlgorithmGlobal.jl")
include("States_NewRungeKutta.jl")
include("NewMyRungeKutta.jl")

#Entrepreneurs Problem
include("ControlsAlgorithmE.jl")
include("States_NewRungeKuttaE.jl")
include("NewMyRungeKuttaE.jl")

include("ProblemFunction.jl")


#Marginal taxes plot
include("marginal_taxes.jl")


#Define values for the model parameters
pa = init_parameters();

    #1. Rawlsian problem:
    #Define initial state vector
    ###ESTEEEEE New Controls+ Bien Ralwsian Planner (Initial boundary condition for μ)

    uw0    = 25 #guess
    μ0     = -1
    e0     = pa.θ_e_lb;
    ϕ_e0   = 1 #guess
    ϕ_e0   = 0.01 #guess
    #ϕ_e0   =  -1.233e-5 #Min- e=1.07
    y_agg0 = 0.0
    λ0     = 1 #guess
    l_agg0 = 0.0
    ω0     = 2 #guess=#
    l_new0 = 0.0
    y_new0 = 0.0
    ystar0 = 0.0
    lstar0 = 0.0

    #Cierra ϕ_e pero no coinciden los estados (Hamiltoniano vs Residuales)
    uw0    = 10.0 #guess
    μ0     = -1.0
    e0     = pa.θ_e_lb;
    ϕ_e0   = 1.0 #guess
    #ϕ_e0   =  -1.233e-5 #Min- e=1.07
    y_agg0 = 0.0
    λ0     = 1.0 #guess
    l_agg0 = 0.0
    ω0     = 3.0 #guess=#
    l_new0 = 0.0
    y_new0 = 0.0
    ystar0 = 0.0
    lstar0 = 0.0

    #Cierra L pero no   ϕ_e
    uw0    = 100.0 #guess
    μ0     = -1.0
    e0     = pa.θ_e_lb;
    ϕ_e0   = 12.0 #guess
    #ϕ_e0   =  -1.233e-5 #Min- e=1.07
    y_agg0 = 0.0
    λ0     = 1.0 #guess
    l_agg0 = 0.0
    ω0     = 3.365 #guess=#
    l_new0 = 0.0
    y_new0 = 0.0
    ystar0 = 0.0
    lstar0 = 0.0

    uw0    = 100.32552083333333#guess
    μ0     = -1.0
    e0     = pa.θ_e_lb;
    ϕ_e0   = 10.026041666666666 #guess
    #ϕ_e0   =  -1.233e-5 #Min- e=1.07
    y_agg0 = 0.0
    λ0     = 1.0 #guess
    l_agg0 = 0.0
    ω0     = 3.37157226562500 #guess=#
    l_new0 = 0.0
    y_new0 = 0.0
    ystar0 = 0.0
    lstar0 = 0.0

    #2. Utilitarian problem:
    #Define initial state vector
    uw0  = 2.983699982905917e-5
    μ0   = 0
    ϕ_e0 = 4.067658551717214
    λ0   = 0.12
    ω0   = 0.12025
    e0     =  pa.θ_e_lb;
    y_agg0 =  0.0
    l_agg0 =  0.0
    l_new0 =  0.0
    y_new0 =  0.0
    ystar0 = 0.0
    lstar0 = 0.0 #n=#

        uw0  = 2.983699982905917e-5
        μ0   = 0
        ϕ_e0 = 4.067658551717214
        λ0   = 0.12
        ω0   = 0.12025
        e0     =  pa.θ_e_lb;
        sse = StateE(uw0, μ0, λ0, ω0)
        θe = e0
        θ = pa.θ_w_lb

#Test my runge kutta
    #Global Problem
        Nspan = 500
        y0= [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0, l_new0, y_new0, 0, 0 ];
        xlb= pa.θ_w_lb;
        xub = pa.θ_w_ub;
        xstep = (xub - xlb)/(Nspan - 1);
        xspan = xlb:xstep:xub;
        solution = Array{Float64}(undef,Nspan,12);
        my_runge_kutta!(solution,y0,xspan,xstep,pa)

        using DelimitedFiles
        writedlm("SolutionNew.csv",solution,';')

        #For the log-normal distribution:
        #=xlb= log(pa.θ_w_lb);
        xub = log(pa.θ_w_ub);
        xstep = (xub - xlb)/(Nspan - 1);
        xspan = xlb:xstep:xub;
        θspan = exp.(xspan);=#

        θspan = Array{Float64,1}
        θspan = collect(xspan)

        controls = Array{Float64}(undef,Nspan,4);
        recover_controls!(controls, θspan, solution);
        controls

        C= DataFrame(controls)
        names!(C,[:z, :n, :l, :p] )
        RungeKutta=hcat(DataFrame(solution),C, DataFrame(theta=θspan))
        CSV.write("RungeKuttaNew.csv", RungeKutta)

        #Final boundary condition
        he_ub= pa.he(pa.θ_w_ub,solution[end,3])
        bound_e= -(pa.indicator*solution[end,1]^pa.ϕ*he_ub+solution[end,6]*(solution[end,3]*controls[end,2]^pa.α -(pa.β/(1.0+pa.σ))*controls[end,1]^(1.0+pa.σ)
                - solution[end,1])*he_ub - solution[end,8]*controls[end,2]*he_ub + solution[end,2]*controls[end,2]^pa.α*(1.0-pa.β*controls[end,1]^pa.σ))
        dist_e= bound_e - solution[end,4]

        #Define Hamiltonian:
        Hamiltonian = Array{Float64}(undef,Nspan,1);

        for i = 1:Nspan
            h_w = pa.hw(θspan[i],solution[i,3]);
            h_e = pa.he(θspan[i],solution[i,3]);
            term1 = pa.indicator*solution[i,1]^pa.ϕ*(h_w+h_e*controls[i,4]);
            term2 = solution[i,2]*pa.χ/θspan[i]*controls[i,3]^(1.0+pa.ψ);
            term3 = solution[i,6]*(solution[i,3]*controls[i,2]^pa.α*controls[i,4]*h_e-pa.β/(1.0+pa.σ)*controls[i,1]^(1.0+pa.σ)*controls[i,4]*h_e-
                    solution[i,1]*(h_w+h_e*controls[i,4])-pa.χ/(1.0+pa.ψ)*controls[i,3]^(1.0+pa.ψ)*h_w);
            term4 = solution[i,8]*(θspan[i]*controls[i,3]*h_w-controls[i,2]*controls[i,4]*h_e)+solution[i,4]*controls[i,4];
            Hamiltonian[i,1] = term1 + term2 + term3 + term4;
        end
        plot(θspan, Hamiltonian[:,1])

    #Entrepreneurs Problem

        #Initial boundary conditions (states from the global problem)
        ue0    =   solution[end,1]
        μe0    =   solution[end,2]
        ye0    =   solution[end,5]
        λe0    =   solution[end,6]
        le0    =   solution[end,7]
        ωe0    =   solution[end,8]

        Nspan = 500
        y0= [ue0, μe0, ye0, λe0, le0, ωe0, 0.0, 0.0];
        elb= solution[end,3];
        eub = pa.θ_e_ub;
        estep = (eub - elb)/(Nspan - 1);
        espan = elb:estep:eub;
        solutione = Array{Float64}(undef,Nspan,10);
        my_runge_kuttae!(solutione,y0,espan,estep,pa,pa.θ_w_ub)

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
        θespan = collect(espan)

        controlse = Array{Float64}(undef,Nspan,2);
        recover_controlse!(controlse, pa.θ_w_ub ,θespan, solutione);
        controlse

        E= DataFrame(controlse)
        names!(E,[:z, :n] )
        CSV.write("ControlsEntrepreneurs.csv", E)
        RungeKuttaE=hcat(DataFrame(solutione),E, DataFrame(thetae=θespan))
        CSV.write("RungeKuttaNewE.csv", RungeKuttaE)

        #Define Hamiltonian:
        HamiltonianE = Array{Float64}(undef,Nspan,1);

        for i = 1:Nspan
            h_e = pa.he(θspan[end],θespan[i]);
            term1 = pa.indicator*solutione[i,1]^pa.ϕ*h_e;
            term2 = solutione[i,4]*(θespan[i]*controlse[i,2]^pa.α-pa.β/(1.0+pa.σ)*controlse[i,1]^(1.0+pa.σ)-
                    solutione[i,1])*h_e;
            term3 = solutione[i,6]*controlse[i,2]*h_e;
            term4 = solutione[i,2]*controlse[i,2]^pa.α*(1.0-pa.β*controlse[i,1]^pa.σ);
            HamiltonianE[i,1] = term1 + term2 - term3 + term4;
        end

        plot(θespan, HamiltonianE[:,1])

#Complete Problem
    uw0    = 100.0 #guess
    uw0    = 1.0 #guess
    #uw0    = 88.1143 #guess
    μ0     = -1.0
    e0     = pa.θ_e_lb;
    ϕ_e0   = 12.0 #guess
    ϕ_e0   = 0.1143
    ϕ_e0   = -0.0007
    y_agg0 = 0.0
    λ0     = 1.0 #guess
    l_agg0 = 0.0
    ω0     = 3.365 #guess=#
    l_new0 = 0.0
    y_new0 = 0.0
    ystar0 = 0.0
    lstar0 = 0.0

    Nspan = 500
    y0= [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0, l_new0, y_new0, 0, 0 ];
    xlb= pa.θ_w_lb;
    xub = pa.θ_w_ub;
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xlb:xstep:xub;
    solution = Array{Float64}(undef,Nspan,12);
    controls = Array{Float64}(undef,Nspan,4);
    solutione = Array{Float64}(undef,Nspan,10);
    controlse = Array{Float64}(undef,Nspan,2);
    complete_problem!(solution,solutione,controls,controlse,y0,xspan,xstep,pa)

    θspan = Array{Float64,1}
    θspan = collect(xspan)
    elb= solution[end,3];
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = elb:estep:eub;
    θespan = Array{Float64,1}
    θespan = collect(espan)

    he_ub= pa.he(pa.θ_w_ub,solution[end,3])
    bound_e= -(pa.indicator*solution[end,1]^pa.ϕ*he_ub+solution[end,6]*(solution[end,3]*controls[end,2]^pa.α -(pa.β/(1.0+pa.σ))*controls[end,1]^(1.0+pa.σ)
            - solution[end,1])*he_ub - solution[end,8]*controls[end,2]*he_ub + solution[end,2]*controls[end,2]^pa.α*(1.0-pa.β*controls[end,1]^pa.σ))
    dist_e= bound_e - solution[end,4]

    #Marginal taxes:
    τ_prime = Array{Float64}(undef,Nspan,3);
    marginal_taxes(controls, θspan, solution, pa)

    τ_prime_e = Array{Float64}(undef,Nspan,2);
    marginal_taxese(controlse, θespan, solutione, pa)

    #Plots:
    graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub, bound_e,τ_prime,τ_prime_e)


    #Uniform distribution
        θ_e_a= pa.μ_e-((12.0^0.5)/2)*(pa.σ2_e^0.5);
        θ_e_b= pa.μ_e+((12.0^0.5)/2)*(pa.σ2_e^0.5);

        #theta_w_a
        θ_w_a= pa.μ_w-((12.0^0.5)/2)*(pa.σ2_w^0.5);
        θ_w_b= pa.μ_w+((12.0^0.5)/2)*(pa.σ2_w^0.5);

        dist_marginal_w=Uniform(θ_w_a,θ_w_b);
        dist_marginal_e=Uniform(θ_e_a,θ_e_b);

dist_marginal_w*dist_marginal_e

#GRAPHS

    pygui(true)
    rc("font", family="serif")

    #Global problem

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

    #Entrepreneurs problem:

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
