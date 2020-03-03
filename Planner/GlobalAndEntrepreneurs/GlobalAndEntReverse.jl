
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
include("NewMyRungeKuttaReverse.jl")

#Entrepreneurs Problem
include("NewControlsAlgorithmE.jl")
include("States_NewRungeKuttaE.jl")
include("NewMyRungeKuttaEReverse.jl")

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

        #Marginal taxes and taxes path:
        τ_prime = Array{Float64}(undef,Nspan,3);
        marginal_taxes(controls, θspan, solution, pa)

        (taxes, taxes_ent) = taxes_path(controls, θspan, solution, pa, θespan, solutione, controlse, τ_prime, τ_prime_e);

        taxes_rev     = fill(NaN,Nspan,6);
        taxes_rev_ent = fill(NaN,Nspan,4);

        taxes_rev[:,1:3] = taxes;
        for j=1:Nspan
            θ = θspan[j];
            e     = solution[j,3];
            ω     = solution[j,8];

            zz    = controls[j,1];
            nn    = controls[j,2];
            ll    = controls[j,3];
            pp    = controls[j,4];

            taxes_rev[j,4] = e*nn^pa.α-ω*nn-taxes[j,2]; #Tc
            taxes_rev[j,5] = ω*nn; #Tn
            taxes_rev[j,6] = θ*ll*ω; #Tl
        end

        taxes_rev_ent[:,1:2] = taxes_ent;
        for j=1:Nspan
            θe = θespan[j];
            ωe    = solutione[j,6];

            nn    = controlse[j,2];

            taxes_rev_ent[j,3] = θe*nn^pa.α-ωe*nn-taxes_ent[j,2]; #Tc
            taxes_rev_ent[j,4] = ωe*nn; #Tn
        end

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

        #Max evasion:
        mat_for_z = Array{Float64}(undef,Nspan,4);

        for i=1:Nspan
            mat_for_z[i,1] = (solution[i,6]*solution[i,3]*controls[i,2]^pa.α-solution[i,8]*controls[i,2]); #Max posible evasion.
            mat_for_z[i,2] = controls[i,1]; #Evasion in model.
            mat_for_z[i,3] = (solutione[i,4]*θespan[i]*controlse[i,2]^pa.α-solutione[i,6]*controlse[i,2]); #Max posible evasion.
            mat_for_z[i,4] = controlse[i,1]; #Evasion in model.
        end

        #Getting the propositions:
        proposition1 = fill(NaN,Nspan,4);
        proposition2 = fill(NaN,Nspan,2);
        proposition3 = fill(NaN,Nspan,6);

        (proposition1, proposition2, proposition3) = propositions(controls,θspan,solution,τ_prime,pa,θespan,solutione);

        #Plots:
        #graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub,
        #       bound_e,τ_prime,τ_prime_e,"C:\\Users\\marya\\Documents\\GitHub\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs\\Graphs",
        #       utilities_prime,A_matrix)
        graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub,
                bound_e,τ_prime,τ_prime_e,"C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs\\Graphs",
                utilities_prime,A_matrix,mat_for_z, proposition1, proposition2, proposition3, taxes_rev, taxes_rev_ent)
