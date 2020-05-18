#cd("C:\\Users\\marya\\Documents\\GitHub\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs")
#cd("C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs")

fig_graphs   = true; #Indicator to print figures (true is when we print figures).
RK_algorithm = true; #Indicator to use RK package from Julia (true is when we use package).

using Roots
using NLopt
using Statistics
fig_graphs && using PyPlot
using DataFrames
using CSV
using NLsolve
RK_algorithm && using DifferentialEquations
import ForwardDiff
#using Plots

#Global Problem
include("Definitions2.jl")
include("MJControlsAlgorithmGlobal3.jl")
RK_algorithm == true ? include("States_NewRungeKutta.jl") : include("original_States.jl")
RK_algorithm == true ? include("NewMyRungeKuttaReverse.jl") : include("original_RK.jl")

#Entrepreneurs Problem
#include("NewControlsAlgorithmE.jl")
include("ControlsAlgorithmEV2.jl")
RK_algorithm == true ? include("States_NewRungeKuttaE.jl") : include("original_StatesEnt.jl")
RK_algorithm == true ? include("NewMyRungeKuttaEReverse.jl") : include("original_RKEnt.jl")

#The file for the function that computes everything:
include("ProblemFunction.jl")

#Marginal taxes:
include("marginal_taxes.jl")

#Propositions:
include("Propositions.jl")
include("Integrals.jl")

#Define values for the model parameters:
pa  = init_parameters();
# Algorithm to solve differencial equations:
#alg = Tsit5()
alg = Rosenbrock23()
#Entrepreneurs Problem
    #Initial boundary conditions (states from the global problem)

    #Define proportion of agents in global problem
    if pa.uniform_ind == true
        #For uniform distribution:
        elb    =   pa.θ_e_ub*(1.0-0.007)
        ue0    =   640.0
        μe0    =   0.0 - 1.0e-10
        ye0    =   0.0
        λe0    =   1.0
        le0    =   0.0
        ωe0    =   1.3 #mbar=1.0   ω_min=1.342555  ω_max=1.342625
        #mbar=10.0: ω_min=1.342555  ω_max=1.339170  -- ω_max_max=1.343049 e(θw_lb)=θe_lb (but bunching not addressed)

        #Utilitarian case:
        elb    =   pa.θ_e_ub*(1.0-0.011)
        ue0    =   340.0
        μe0    =   0.0 - 1.0e-10
        ye0    =   0.0
        λe0    =   0.04
        le0    =   0.0
        ωe0    =   0.039 #mbar=1.0   ω_min=1.342555  ω_max=1.342625

        #Utilitarian case:
        elb    =   pa.θ_e_ub*(1.0-0.011)
        ue0    =   462.0
        μe0    =   0.0 - 1.0e-10
        ye0    =   0.0
        λe0    =   0.07
        le0    =   0.0
        ωe0    =   0.0929 #mbar=1.0   ω_min=1.342555  ω_max=1.342625

    else
        #For log-normal distribution:
        ue0 =   3000.0
        μe0 =   0.0 - 1.0e-10
        ye0 =   0.0
        λe0 =   1.0
        le0 =   0.0
        ωe0 =   1.2
        elb =   pa.θ_e_ub*(1.0-0.2)
        #mbar=10.0: ω_min=1.342555  ω_max=1.339170  -- ω_max_max=1.343049 e(θw_lb)=θe_lb (but bunching not addressed)
    end

    #For log-normal distribution:
    #gp     =   1.0-0.23
    #ue0    =   290.0
    #ωe0    =   0.76999
    #mbar=10.0: ω_min=1.342555  ω_max=1.339170  -- ω_max_max=1.343049 e(θw_lb)=θe_lb (but bunching not addressed)

    Nespan = 500
    y_end= [ue0, μe0, ye0, λe0, le0, ωe0, 0.0, 0.0];
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = collect(elb:estep:eub);
    statese = Array{Float64}(undef,Nspan,10);
    fill!(statese,NaN);

    if RK_algorithm == true
        my_runge_kuttae_reverse!(statese,espan,y_end,estep,pa,alg)
    else
        my_runge_kuttae_reverse!(statese,y_end,espan,estep,pa,pa.θ_w_ub)
    end

    #statese[end,:]
    #using DelimitedFiles
    #writedlm("SolutionNewE.csv",states,';')

    #For log-normal distribution
    #=xlb= log(pa.θ_w_lb);
    xub = log(pa.θ_w_ub);
    xstep = (xub - xlb)/(Nespan - 1);
    xspan = xlb:xstep:xub;
    θspan = exp.(xspan);=#

    θespan = Array{Float64,1};
    θespan = collect(elb:estep:eub);

    controlse = Array{Float64}(undef,Nespan,2);
    fill!(controlse,NaN);
    recover_controlse!(controlse, θespan, statese);
    controlse

    #E= DataFrame(controlse)
    #names!(E,[:z, :n] )
    #CSV.write("ControlsEntrepreneurs.csv", E)
    #RungeKuttaE=hcat(DataFrame(statese),E, DataFrame(thetae=θespan))
    #CSV.write("RungeKuttaNewE.csv", RungeKuttaE)

    #Marginal taxes and taxes path:
    τ_prime_e = Array{Float64}(undef,Nespan,2);
    fill!(τ_prime_e,NaN);
    marginal_taxese(controlse, θespan, statese, pa);

    #taxese = DataFrame(τ_prime_e)
    #names!(taxese,[:tau_c, :tau_n] )
    #taxes2e=hcat(DataFrame(theta=θespan),taxese)
    #CSV.write("marginal_taxes_e.csv",taxes2e)

#Global Problem (Reverse)
    uw_end    = statese[1,1] #guess
    μ_end     = statese[1,2]
    e_end     = elb;
    he_ub     = pa.he(pa.θ_w_ub,e_end)
    ϕ_e_end   = -(pa.indicator*statese[1,1]^pa.ϕ*he_ub+statese[1,4]*(e_end*controlse[1,2]^pa.α -(pa.β/(1.0+pa.σ))*controlse[1,1]^(1.0+pa.σ)
                - statese[1,1])*he_ub - statese[1,6]*controlse[1,2]*he_ub + statese[1,2]*controlse[1,2]^pa.α*(1.0-pa.β*controlse[1,1]^pa.σ))
    y_agg_end = statese[1,3]
    λ_end     = 1.0 #guess
    l_agg_end = statese[1,5]
    ω_end     = statese[1,6] #guess=#
    l_new_end = 0.0
    y_new_end = 0.0
    ystar_end = 0.0
    lstar_end = 0.0

    Nspan = 500
    y_end= [uw_end, μ_end, e_end, ϕ_e_end, y_agg_end, λ_end, l_agg_end, ω_end, l_new_end, y_new_end, 0, 0 ];
    xlb = pa.θ_w_lb+0.001;
    xub = pa.θ_w_ub;
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = collect(xlb:xstep:xub)
    states = Array{Float64}(undef,Nspan,12);
    fill!(states,NaN)
    if RK_algorithm == true
        lenght_sol=my_runge_kutta_reverse!(states, xspan, y_end,xstep,pa,alg)
    else
        my_runge_kutta_reverse!(states,y_end,xspan,xstep,pa)
        lenght_sol=Nspan
    end
    #display(xspan)
    #display(pa.θ_e_lb)
    #using DelimitedFiles
    #writedlm("SolutionNew.csv",states,';')

    θspan = Array{Float64,1};
    θspan = xspan;

    controls = Array{Float64}(undef,Nspan,4);

    recover_controls!(controls, θspan, states);

    bound_e   = -(pa.indicator*statese[1,1]^pa.ϕ*he_ub+statese[1,4]*(e_end*controlse[1,2]^pa.α -(pa.β/(1.0+pa.σ))*controlse[1,1]^(1.0+pa.σ)
                - statese[1,1])*he_ub - statese[1,6]*controlse[1,2]*he_ub + statese[1,2]*controlse[1,2]^pa.α*(1.0-pa.β*controlse[1,1]^pa.σ))

    #C= DataFrame(controls)
    #names!(C,[:z, :n, :l, :p] )
    #RungeKutta=hcat(DataFrame(states),C, DataFrame(theta=θspan))
    #CSV.write("RungeKuttaNew.csv", RungeKutta)

        #Marginal taxes and taxes path:
        τ_prime = Array{Float64}(undef,Nspan,3);
        fill!(τ_prime,NaN);
        marginal_taxes(controls, θspan, states, pa)

        (taxes, taxes_ent, taxes_rev, taxes_rev_ent) = taxes_path(controls, θspan, states, pa, θespan, statese, controlse, τ_prime, τ_prime_e);

        taxes_full     = fill(NaN,Nspan,6);
        taxes_ent_full = fill(NaN,Nspan,4);
        taxes_full[:,1:3] .= taxes;
        taxes_full[:,4:6] .= taxes_rev;
        taxes_ent_full[:,1:2] .= taxes_ent;
        taxes_ent_full[:,3:4] .= taxes_rev_ent;

        #taxes = DataFrame(τ_prime)
        #names!(taxes,[:tau_c, :tau_n, :tau_l] )
        #taxes2=hcat(DataFrame(theta=θspan),taxes)
        #CSV.write("marginal_taxes.csv",taxes2)

        #Utilities:
        utilities_prime = Array{Float64}(undef,Nspan,2);
        fill!(utilities_prime,NaN);

        for i=1:Nspan
            utilities_prime[i,1] = pa.χ*controls[i,3]^(1.0+pa.ψ)/θspan[i]; #u_w
            utilities_prime[i,2] = controls[i,2]^pa.α*(1.0-pa.β*controls[i,1]^pa.σ); #u_e
        end

        #A:
        A_matrix = Array{Float64}(undef,Nspan,2);
        fill!(A_matrix,NaN);

        for i=1:Nspan
            A_matrix[i,1] = states[i,8]*pa.ς + pa.indicator*states[i,1]^pa.ϕ-states[i,6]*states[i,1] + states[i,4]/pa.he(θspan[i],states[i,3]); #A
            n_full_info   = ((states[i,6]*pa.α*states[i,3])/states[i,8])^(1.0/(1.0-pa.α));
            A_matrix[i,2] = -(1.0-pa.α)/pa.α*states[i,8]*n_full_info; #A if z is 0 and n_full_info
        end

        #Max evasion:
        mat_for_z = Array{Float64}(undef,Nspan,4);
        fill!(mat_for_z,NaN);

        for i=1:Nspan
            mat_for_z[i,1] = states[i,6]*states[i,3]*controls[i,2]^pa.α-states[i,8]*(controls[i,2]-pa.ς); #Max posible evasion.
            mat_for_z[i,2] = controls[i,1]; #Evasion in model.
            mat_for_z[i,3] = statese[i,4]*θespan[i]*controlse[i,2]^pa.α-statese[i,6]*(controlse[i,2]-pa.ς); #Max posible evasion.
            mat_for_z[i,4] = controlse[i,1]; #Evasion in model.
        end

        #Getting the propositions:
        proposition1 = fill(NaN,Nspan,4);
        proposition2 = fill(NaN,Nspan,2);
        proposition3 = fill(NaN,Nspan,6);

        (proposition1, proposition2, proposition3) = propositions(controls, θspan, states, τ_prime, pa, θespan, statese);

#ctrlvec = controls
#θvec = θspan
#solvec = states
#controls,θspan,states,τ_prime,pa,θespan,statese

        #Plots:
        fig_graphs && graphs!(states,statese,controls,controlse, θspan, θespan, pa.θ_w_ub,
                            bound_e,τ_prime,τ_prime_e,".\\Graphs",
                            utilities_prime,A_matrix,mat_for_z, proposition1, proposition2, proposition3, taxes_full, taxes_ent_full)



                diff = fill(NaN,Nspan,2)
                for j = 1:Nspan
                    diff[j,1] = proposition3[j,1]+ proposition3[j,2]-proposition3[j,3]
                    #diff[j,2] = proposition1[j,1]- proposition1[j,2]-proposition1[j,3]
                    diff[j,2] = -states[j,2]- proposition1[j,2]-proposition1[j,3]
                end

    if fig_graphs

                        #Proposition3:
                        fig, prop3=plt.subplots(1,2)
                        fig.suptitle("Proposition 3")
                            prop3[1].plot(θspan[:], proposition3[:,1])
                            prop3[1].plot(θspan[:], proposition3[:,2])
                            prop3[1].plot(θspan[:], proposition3[:,3])
                            prop3[1].legend(["εz z/n^α he","1/λ [Ve-Vw] g 1/ue'","λ he p"],loc="upper right")
                                prop3[2].plot(θspan[:], diff[:,1])

                                #Proposition1:
                                fig, prop1=plt.subplots(1,2)
                                fig.suptitle("Proposition 1")
                                    #prop1[1].plot(θspan[:], proposition1[:,1])
                                    prop1[1].plot(θspan[:], -states[:,2])
                                    prop1[1].plot(θspan[:], proposition1[:,2])
                                    prop1[1].plot(θspan[:], proposition1[:,3])
                                    prop1[1].legend(["-μ","Tl'/(1-Tl') εl/(1+εl) θw hw","εz z/n^α he"],loc="upper right")
                                    prop1[2].plot(θspan[:], diff[:,2])

                                    # diff()

        to_plot = fill(NaN,Nspan,3);
        for j=1:Nspan
            #Definition of states and controls:
            θ = θspan[j];

                  uw    = states[j,1];
                  μ     = states[j,2];
                  e     = states[j,3];
                  ϕ_e   = states[j,4];
                  y_agg = states[j,5];
                  λ     = states[j,6];
                  l_agg = states[j,7];
                  ω     = states[j,8];

                  zz    = controls[j,1];
                  nn    = controls[j,2];
                  ll    = controls[j,3];
                  pp    = controls[j,4];

                  Tn = taxes_full[j,2];

                  #Find the marginal taxes:
                  to_plot[j,1] = e*nn^pa.α-ω*(nn-pa.ς)
                  to_plot[j,2] = e*nn^pa.α-ω*(nn-pa.ς) - Tn
                  to_plot[j,3] = e*nn^pa.α-ω*(nn-pa.ς) - Tn - zz

                end


                        fig, tax_ent=plt.subplots(2,3)
                        fig.suptitle("")
                            #τ_c:
                        tax_ent[1,1].plot(θspan[1:500], to_plot[:,1])
                        tax_ent[1,1].set(ylabel="e*nn^pa.α-ω*e*nn", xlabel="θw")
                            #τ_n:
                        tax_ent[1,2].plot(θspan[1:500], to_plot[:,2])
                        tax_ent[1,2].set(ylabel="e*nn^pa.α-ω*e*nn - Tn", xlabel="θw")
                            #τ_c:
                        tax_ent[1,3].plot(θspan[1:500], to_plot[:,3])
                        tax_ent[1,3].set(ylabel="e*nn^pa.α-ω*e*nn - Tn - z", xlabel="θw")
                        #τ_c:
                    tax_ent[2,1].plot(states[:,3], to_plot[:,1])
                    tax_ent[2,1].set(ylabel="e*nn^pa.α-ω*e*nn", xlabel="e")
                        #τ_n:
                    tax_ent[2,2].plot(states[:,3], to_plot[:,2])
                    tax_ent[2,2].set(ylabel="e*nn^pa.α-ω*e*nn - Tn", xlabel="e")
                        #τ_c:
                    tax_ent[2,3].plot(states[:,3], to_plot[:,3])
                    tax_ent[2,3].set(ylabel="e*nn^pa.α-ω*e*nn - Tn - z", xlabel="e")

    end

                    # graphs = fill(NaN,Nspan,2)
                    # for j = 1:Nspan
                    #     graphs[j,1] = controls[j,2]^pa.α
                    # end

                    # graphs[2:500,2] = diff(controls[:,1])./diff(states[:,3])

                    #                         fig, tax_ent=plt.subplots(1,2)
                    #                         fig.suptitle("")
                    #                             #τ_c:
                    #                         tax_ent[1].plot(states[:,3], graphs[:,1],states[:,3], graphs[:,2])
                    #                             #τ_n:
                    #                         tax_ent[2].plot(θspan[1:500], graphs[:,2])



                    #                             #τ_c:
                    #                         tax_ent[1,3].plot(θspan[1:500], to_plot[:,3])
                    #                         tax_ent[1,3].set(ylabel="e*nn^pa.α-ω*e*nn - Tn - z", xlabel="θw")
                    #                         #τ_c:
                    #                     tax_ent[2,1].plot(states[:,3], to_plot[:,1])
                    #                     tax_ent[2,1].set(ylabel="e*nn^pa.α-ω*e*nn", xlabel="e")
                    #                         #τ_n:
                    #                     tax_ent[2,2].plot(states[:,3], to_plot[:,2])
                    #                     tax_ent[2,2].set(ylabel="e*nn^pa.α-ω*e*nn - Tn", xlabel="e")
                    #                         #τ_c:
                    #                     tax_ent[2,3].plot(states[:,3], to_plot[:,3])
                    #                     tax_ent[2,3].set(ylabel="e*nn^pa.α-ω*e*nn - Tn - z", xlabel="e")
