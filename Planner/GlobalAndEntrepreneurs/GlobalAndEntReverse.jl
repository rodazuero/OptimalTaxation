
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

#Marginal taxes:
include("marginal_taxes.jl")

#Propositions:
include("Propositions.jl")

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
    ωe0    =   1.34297

    Nspan = 500
    y_end= [ue0, μe0, ye0, λe0, le0, ωe0, 0.0, 0.0];
    elb = pa.θ_e_ub - ((1-gp)*(pa.θ_e_ub-pa.θ_e_a)*(1.0-pa.constant_w_lw*pa.constant_e_lw));
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

    #Marginal taxes and taxes path:
    τ_prime_e = Array{Float64}(undef,Nspan,2);
    marginal_taxese(controlse, θespan, solutione, pa);

    TaxesE = Array{Float64}(undef,Nspan,2);
    taxes_pathe(controlse, θespan, solutione, pa)

    fig, tax_e=plt.subplots(1,2)
    fig.suptitle("Taxes Path")
        #τ_c:
    tax_e[1].plot(θespan[1:500], TaxesE[:,1])
    tax_e[1].set(ylabel="τ_c")
        #τ_n:
    tax_e[2].plot(θespan[1:500], TaxesE[:,2])
    tax_e[2].set(ylabel="τ_n")


    (θespan[2]-θespan[1])/2.0*(pa.β*controlse[1,1]^pa.σ + pa.β*controlse[2,1]^pa.σ)
    #taxese = DataFrame(τ_prime_e)
    #names!(taxese,[:tau_c, :tau_n] )
    #taxes2e=hcat(DataFrame(theta=θespan),taxese)
    #CSV.write("marginal_taxes_e.csv",taxes2e)

    taxese = DataFrame(TaxesE)
    names!(TaxesE,[:tau_c, :tau_n] )
    taxes2e=hcat(DataFrame(theta=θespan),taxese)
    CSV.write("taxes_e.csv",taxes2e)

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

        #Marginal taxes and taxes path:
        τ_prime = Array{Float64}(undef,Nspan,3);
        marginal_taxes(controls, θspan, solution, pa)

        taxes_path(controls, θspan, solution, pa);

        ctrlvec = controls
        θvec = θspan
        solvec = solution


        fig, tax=plt.subplots(1,3)
        fig.suptitle("Taxes")
            #τ_c:
        tax[1,1].plot(θspan[1:500], Taxes[:,1])
        tax[1,1].set(ylabel="τ_c'")
            #τ_n:
        tax[2].plot(θspan[1:500], Taxes[:,2])
        tax[2].set(ylabel="τ_n'")
            #τ_l:
        tax[3].plot(θspan[1:500], Taxes[:,3])
        tax[3].set(ylabel="τ_l'")

θ*ω-λ*pa.χ*ll^pa.ψ
θspan[1]*solution[1,8]-pa.χ*controls[1,3]^pa.ψ
θspan[1]*solution[1,8]/solution[1,6]*controls[1,3]-pa.χ/(1.0+pa.ψ)*controls[1,3]^(1.0+pa.ψ)-solution[1,1]
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

                proposition1 = Array{Float64}(undef,Nspan,3)
                proposition2 = Array{Float64}(undef,Nspan,2)
                proposition3 = Array{Float64}(undef,Nspan,4)
        (proposition1, proposition2, proposition3) = propositions(controls,θspan,solution,τ_prime,pa,θespan,solutione)

int = proposition1[:,1]
(sum(int)-0.5*int[1]-0.5*int[end])*(θspan[end]-θspan[1])/(Nspan-1.0)

*(1.0-pa.constant_w_lw*pa.constant_e_lw);

int1 = proposition3[:,3]
int2 = proposition3[:,4]

sum(int1)
sum(int2)
(sum(int2)-0.5*int2[1]-0.5*int2[end])*(θspan[end]-θspan[1])/(Nspan-1.0)*(1.0-pa.constant_w_lw*pa.constant_e_lw);


mus = Array{Float64}(undef,Nspan,2)

for i =1:Nspan
    θ = θspan[i]
    e = solution[i,3]

    mus[i,1] = solution[i,2]
    mus[i,2] = (1.0-((e-pa.θ_e_a)*(θ-pa.θ_w_a))/((pa.θ_e_ub-pa.θ_e_a)*(pa.θ_w_ub-pa.θ_w_a)))/0.9999 + solution[i,2]
end

plot(mus[:,2])
