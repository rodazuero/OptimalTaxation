
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
    gp     =   0.9

    ue0    =   1000.0
    ue0    =   900.0
    μe0    =   0.0 - 1.0e-10
    ye0    =   0.0
    λe0    =   1.0
    le0    =   0.0
    ωe0    =   0.9
    ωe0    =   1.272054
    ωe0    =   1.272

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

    #Marginal taxes:
    τ_prime_e = Array{Float64}(undef,Nspan,2);

    marginal_taxese(controlse, θespan, solutione, pa)

    taxese = DataFrame(τ_prime_e)
    names!(taxese,[:tau_c, :tau_n] )
    taxes2e=hcat(DataFrame(theta=θespan),taxese)
    CSV.write("marginal_taxes_e.csv",taxes2e)

    #Plots:
    #Entrepreneurs problem:

    #States
    fig, estados_e=plt.subplots(3,2)
    fig.suptitle("Optimal States - Entrepreneurs Problem")
        #u_e:
    estados_e[1,1].plot(θespan[1:500], solutione[1:500,1])
    #estados_e[1,1].set_title("u_e")
    estados_e[1,1].set(ylabel="u_e")
        #μ_e:
    estados_e[1,2].plot(θespan[1:500], solutione[1:500,2])
    estados_e[1,2].plot(θespan[1:500], repeat([0],500), "tab:green")
    #estados_e[1,2].set_title("μe")
    estados_e[1,2].set(ylabel="μe")
        #Y_e:
    estados_e[2,1].plot(θespan[1:500], solutione[1:500,3])
    estados_e[2,1].plot(θespan[1:500], repeat([0.15],500), "tab:green")
    #estados_e[2,1].set_title("Ye")
    estados_e[2,1].set(ylabel="Ye")
        #λ_e:
    estados_e[2,2].plot(θespan[1:500], solutione[1:500,4])
    #estados_e[2,2].set_title("λe")
    estados_e[2,2].set(ylabel="λe")
        #L_e:
    estados_e[3,1].plot(θespan[1:500], solutione[1:500,5])
    estados_e[3,1].plot(θespan[1:500], repeat([0],500), "tab:green")
    #estados_e[3,1].set_title("Le")
    estados_e[3,1].set(ylabel="Le")
        #ω:
    estados_e[3,2].plot(θespan[1:500], solutione[1:500,6])
    #estados_e[3,2].set_title("ωe")
    estados_e[3,2].set(ylabel="ωe")

    #Auxiliar states
    fig, estados_auxE=plt.subplots(1,2)
    fig.suptitle("Auxiliar States - Global Problem")
        #L*:
    estados_auxE[1,1].plot(θespan[1:500], solutione[1:500,9])
    estados_auxE[1,1].plot(θespan[1:500], repeat([0],500), "tab:green")
    #estados_auxE[1,1].set_title("L*")
    estados_auxE[1,1].set(ylabel="L*")
        #Y*:
    estados_auxE[2].plot(θespan[1:500], solutione[1:500,10])
    estados_auxE[2].plot(θespan[1:500], repeat([0.15],500), "tab:green")
    #estados_auxE[2].set_title("Y*")
    estados_auxE[2].set(ylabel="Y*")

    #Controls
    fig, controles_e=plt.subplots(1,2)
    fig.suptitle("Optimal Controls - Entrepreneurs Problem")
        #ze:
    controles_e[1,1].plot(θespan[1:500], controlse[:,1])
    #controles_e[1,1].set_title("ze")
    controles_e[1,1].set(ylabel="ze")
        #ne:
    controles_e[2].plot(θespan[1:500], controlse[:,2])
    #controles_e[1].set_title("ne")
    controles_e[2].set(ylabel="ne")

    #Marginal taxes:
    fig, margtax_e=plt.subplots(1,2)
    fig.suptitle("Marginal Taxes")
        #τ_c_prime:
    margtax_e[1].plot(θespan[1:500], τ_prime_e[:,1])
    #margtax_e[1].set_title("τ_c'")
    margtax_e[1].set(ylabel="τ_c'")
        #τ_n_prime:
    margtax_e[2].plot(θespan[1:500], τ_prime_e[:,2])
    #margtax_e[2].set_title("τ_n'")
    margtax_e[2].set(ylabel="τ_n'")

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

        #Marginal taxes:
        τ_prime = Array{Float64}(undef,Nspan,3);
        marginal_taxes(controls, θspan, solution, pa)

        taxes = DataFrame(τ_prime)
        names!(taxes,[:tau_c, :tau_n, :tau_l] )
        taxes2=hcat(DataFrame(theta=θspan),taxes)
        CSV.write("marginal_taxes.csv",taxes2)

        #Utilities:
        utilities_prime = Array{Float64}(undef,Nspan,2);

        for i=1:Nspan
            utilities_prime[i,1] = pa.χ*controls[i,3]^(1.0+pa.ψ)/θspan[i]; #u_w
            utilities_prime[i,2] = controlse[i,2]^pa.α*(1.0-pa.β*controlse[i,1]^pa.σ); #u_e
        end

        #A:
        A_matrix = Array{Float64}(undef,Nspan,2);

        for i=1:Nspan
            A_matrix[i,1] = (pa.indicator*solution[i,1]^pa.ϕ-solution[i,6]*solution[i,1]) + solution[i,4]/pa.he(θspan[i],solution[i,3]); #A
            n_full_info   = ((solution[i,6]*pa.α*solution[i,3])/solution[i,8])^(1.0/(1.0-pa.α));
            A_matrix[i,2] = -(1.0-pa.α)/pa.α*solution[i,8]*n_full_info; #A if z is 0 and n_full_info
        end

        plot(solution[:,3], controls[:,2].^pa.α.*(1.0.-τ_prime[:,1]))
        plot(solution[:,3], A_matrix[:,1])

        #Plots:
        #graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub,
        #        bound_e,τ_prime,τ_prime_e,"C:\\Users\\marya\\Documents\\GitHub\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs\\Graphs")
        graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub,
                bound_e,τ_prime,τ_prime_e,"C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs\\Graphs",
                utilities_prime,A_matrix)
