cd("C:\\Users\\marya\\Documents\\GitHub\\OptimalTaxation\\Planner\\InformalityGlobalAndEntrepreneurs")
cd("C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\InformalityGlobalAndEntrepreneurs")

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
include("Definitions.jl")
#include("MJControlsAlgorithmGlobal3.jl")
#include("States_NewRungeKutta.jl")
#include("NewMyRungeKuttaReverse.jl")

#Entrepreneurs Problem
#include("NewControlsAlgorithmE.jl")
include("NewControlsAlgorithmE_Grid.jl")
include("States_NewRungeKuttaE.jl")
include("NewMyRungeKuttaEReverse.jl")
#Entrepreneurs problem without informality
include("NotInfRungeKuttaEntrepreneurs.jl")
include("NotInfStatesEntrepreneurs.jl")
include("NotInfControlsEntrepreneurs.jl")
include("NotInfTaxes.jl")

#The file for the function that computes everything:
include("ProblemFunction.jl")

#Taxes:
include("Taxes.jl")

#Propositions:
#include("Propositions.jl")
#include("Integrals.jl")

#Define values for the model parameters
pa = init_parameters();

#Entrepreneurs Problem
    #Initial boundary conditions (states from the global problem)

    #Define proportion of agents in global problem
    #gp     =   0.993
    gp     =   0.8

    ue0    =   640.0 #guess
    μe0    =   0.0 - 1.0e-10
    ye0    =   0.0
    λe0    =   1.0
    lfe0   =   0.0
    ωfe0   =   1.342 #guess
    lie0   =   0.0
    ωie0   =   1.342 #guess
    wie0   =   1.342 #guess
    ϕwe0   =   0.0

    Nspan = 500
    y_end = [ue0, μe0, ye0, λe0, lfe0, ωfe0, lie0, ωie0, wie0, ϕwe0, 0.0, 0.0, 0.0];
    elb = pa.θ_e_ub - ((1-gp)*(pa.θ_e_ub-pa.θ_e_a)*(1.0-pa.constant_w_lw*pa.constant_e_lw));
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = eub:-estep:elb;
    solutione = Array{Float64}(undef,Nspan,16);
    fill!(solutione,NaN);
    @time my_runge_kuttae_reverse!(solutione,y_end,espan,estep,pa,pa.θ_w_ub)

    θespan = Array{Float64,1}
    θespan = collect(elb:estep:eub)

    controlse = Array{Float64}(undef,Nspan,3);
    fill!(controlse,NaN);
    recover_controlse!(controlse, pa.θ_w_ub ,θespan, solutione);
    controlse

    #States
    fig, estados_e=plt.subplots(5,2)
    fig.suptitle("Optimal States - Entrepreneurs Problem")
        #ue:
    estados_e[1,1].plot(θespan[:], solutione[:,1])
    #estados_e[1,1].set_title("u_e")
    estados_e[1,1].set(ylabel="ue")
        #μe:
    estados_e[1,2].plot(θespan[:], solutione[:,2])
    estados_e[1,2].plot(θespan[:], repeat([0],500), "tab:green")
    #estados_e[1,2].set_title("μe")
    estados_e[1,2].set(ylabel="μe")
        #Ye:
    estados_e[2,1].plot(θespan[:], solutione[:,3])
    estados_e[2,1].plot(θespan[:], repeat([0.15],500), "tab:green")
    estados_e[2,1].set(ylabel="Ye")
        #λe:
    estados_e[2,2].plot(θespan[:], solutione[:,4])
    estados_e[2,2].set(ylabel="λe")
        #Lfe:
    estados_e[3,1].plot(θespan[:], solutione[:,5])
    estados_e[3,1].plot(θespan[:], repeat([0],500), "tab:green")
    #estados_e[3,1].set_title("Le")
    estados_e[3,1].set(ylabel="Lfe")
        #ωf:
    estados_e[3,2].plot(θespan[:], solutione[:,6])
    estados_e[3,2].set(ylabel="ωfe")
        #Li_e:
    estados_e[4,1].plot(θespan[:], solutione[:,7])
    estados_e[4,1].plot(θespan[:], repeat([0],500), "tab:green")
    estados_e[4,1].set(ylabel="Lie")
        #ωie:
    estados_e[4,2].plot(θespan[:], solutione[:,8])
    estados_e[4,2].set(ylabel="ωie")
        #wie:
    estados_e[5,1].plot(θespan[:], solutione[:,9])
    estados_e[5,1].plot(θespan[:], repeat([0],500), "tab:green")
    estados_e[5,1].set(ylabel="wie")
        #ϕwe:
    estados_e[5,2].plot(θespan[:], solutione[:,10])
    estados_e[5,2].set(ylabel="φwe")

    #Controls
    fig, controles_e=plt.subplots(1,3)
    fig.suptitle("Optimal Controls - Entrepreneurs Problem")
        #ze:
    controles_e[1].plot(θespan[:], controlse[:,1])
    controles_e[1].set(ylabel="ze")
        #ne:
    controles_e[2].plot(θespan[:], controlse[:,2])
    controles_e[2].set(ylabel="ne")
    #controles_e[2].legend(["n","n full info"],loc="upper left")
        #nie:
    controles_e[3].plot(θespan[:], controlse[:,3])
    controles_e[3].set_title("nie")

    #Marginal taxes:
    τ_prime_e = Array{Float64}(undef,Nspan,3);
    fill!(τ_prime_e,NaN);
    marginal_taxese!(τ_prime_e,solutione,controlse, θespan, pa);

    #Marginal taxes:
    fig, margtax_e=plt.subplots(1,3)
    fig.suptitle("Marginal Taxes")
        #τc prime:
    margtax_e[1].plot(θespan[:], τ_prime_e[:,1])
    margtax_e[1].set(ylabel="τ_c'")
        #τn prime:
    margtax_e[2].plot(θespan[:], τ_prime_e[:,2])
    margtax_e[2].set(ylabel="τ_n'")
        #τn prime:
    margtax_e[3].plot(θespan[:], τ_prime_e[:,3])
    margtax_e[3].set(ylabel="τ_n' (eq ni)")

    #Solving entrepreneurs without informality:
    Nspan = 500
    y_end= [ue0, μe0, ye0, λe0, lfe0, ωfe0, 0.0, 0.0];
    elb = pa.θ_e_ub - ((1-gp)*(pa.θ_e_ub-pa.θ_e_a)*(1.0-pa.constant_w_lw*pa.constant_e_lw));
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = eub:-estep:elb;
    NotInfsolutione = Array{Float64}(undef,Nspan,10);
    fill!(NotInfsolutione,NaN);
    NotInfmy_runge_kuttae_reverse!(NotInfsolutione,y_end,espan,estep,pa,pa.θ_w_ub)

    θespan = Array{Float64,1}
    θespan = collect(elb:estep:eub)

    NotInfcontrolse = Array{Float64}(undef,Nspan,2);
    fill!(NotInfcontrolse,NaN);
    NotInfrecover_controlse!(NotInfcontrolse, pa.θ_w_ub ,θespan, solutione);
    NotInfcontrolse

    #Gráficos comparando:
    #States
    fig, estados_e=plt.subplots(5,2)
    fig.suptitle("Optimal States - Entrepreneurs Problem")
        #ue:
    estados_e[1,1].plot(θespan[:], solutione[:,1],θespan[:], NotInfsolutione[:,1])
    estados_e[1,1].legend(["ue with informality", "ue without informality"],loc="upper right")
    estados_e[1,1].set(ylabel="ue", xlabel = "θe")
        #μe:
    estados_e[1,2].plot(θespan[:], solutione[:,2], θespan[:], NotInfsolutione[:,2])
    estados_e[1,2].plot(θespan[:], repeat([0],500), "tab:green")
    estados_e[1,2].legend(["μe with informality", "μe without informality"],loc="upper right")
    estados_e[1,2].set(ylabel="μe", xlabel = "θe")
        #Ye:
    estados_e[2,1].plot(θespan[:], solutione[:,3], θespan[:], NotInfsolutione[:,3])
    estados_e[2,1].plot(θespan[:], repeat([0.15],500), "tab:green")
    estados_e[2,1].legend(["Ye with informality", "Ye without informality"],loc="upper right")
    estados_e[2,1].set(ylabel="Ye", xlabel = "θe")
        #λe:
    estados_e[2,2].plot(θespan[:], solutione[:,4], θespan[:], NotInfsolutione[:,4])
    estados_e[2,2].legend(["λe with informality", "λe without informality"],loc="upper right")
    estados_e[2,2].set(ylabel="λe", xlabel = "θe")
        #Lfe:
    estados_e[3,1].plot(θespan[:], solutione[:,5], θespan[:], NotInfsolutione[:,5])
    estados_e[3,1].plot(θespan[:], repeat([0],500), "tab:green")
    estados_e[3,1].legend(["Lfe with informality", "Lfe without informality"],loc="upper right")
    estados_e[3,1].set(ylabel="Lfe", xlabel = "θe")
        #ωf:
    estados_e[3,2].plot(θespan[:], solutione[:,6], θespan[:], NotInfsolutione[:,6])
    estados_e[3,2].legend(["ωfe with informality", "ωfe without informality"],loc="upper right")
    estados_e[3,2].set(ylabel="ωfe", xlabel = "θe")
        #Li_e:
    estados_e[4,1].plot(θespan[:], solutione[:,7])
    estados_e[4,1].plot(θespan[:], repeat([0],500), "tab:green")
    estados_e[4,1].set(ylabel="Lie", xlabel = "θe")
        #ωie:
    estados_e[4,2].plot(θespan[:], solutione[:,8])
    estados_e[4,2].set(ylabel="ωie", xlabel = "θe")
        #wie:
    estados_e[5,1].plot(θespan[:], solutione[:,9])
    estados_e[5,1].plot(θespan[:], repeat([0],500), "tab:green")
    estados_e[5,1].set(ylabel="wie", xlabel = "θe")
        #ϕwe:
    estados_e[5,2].plot(θespan[:], solutione[:,10])
    estados_e[5,2].set(ylabel="φwe", xlabel = "θe")

    #Controls
    fig, controles_e=plt.subplots(1,3)
    fig.suptitle("Optimal Controls - Entrepreneurs Problem")
        #ze:
    controles_e[1].plot(θespan[:], controlse[:,1], θespan[:], NotInfcontrolse[:,1])
    controles_e[1].legend(["ze with informality", "ze without informality"],loc="upper right")
    controles_e[1].set(ylabel="ze", xlabel = "θe")
        #ne:
    controles_e[2].plot(θespan[:], controlse[:,2], θespan[:], NotInfcontrolse[:,2])
    controles_e[1].legend(["ne with informality", "ne without informality"],loc="upper right")
    controles_e[2].set(ylabel="ne", xlabel = "θe")
    #controles_e[2].legend(["n","n full info"],loc="upper left")
        #nie:
    controles_e[3].plot(θespan[:], controlse[:,3])
    controles_e[3].set(ylabel="nie", xlabel = "θe")

    #Marginal taxes:
    NIτ_prime_e = Array{Float64}(undef,Nspan,2);
    fill!(NIτ_prime_e,NaN);
    NotInfmarginal_taxese(NIτ_prime_e,NotInfcontrolse,θespan,NotInfsolutione, pa);

    #Marginal taxes:
    fig, margtax_e=plt.subplots(1,3)
    fig.suptitle("Marginal Taxes")
        #τc prime:
    margtax_e[1].plot(θespan[:], τ_prime_e[:,1],θespan[:], NIτ_prime_e[:,1])
    margtax_e[1].legend(["τ_c' with informality", "τ_c' without informality"],loc="upper right")
    margtax_e[1].set(ylabel="τ_c'",xlabel = "θe")
        #τn prime:
    margtax_e[2].plot(θespan[:], τ_prime_e[:,2],θespan[:], NIτ_prime_e[:,2])
    margtax_e[2].legend(["τ_n' with informality", "τ_n' without informality"],loc="upper right")
    margtax_e[2].set(ylabel="τ_n'", xlabel = "θe")
        #τn prime:
    margtax_e[3].plot(θespan[:], τ_prime_e[:,3])
    margtax_e[3].set(ylabel="τ_n' (eq ni)")

    difference = Array{Float64}(undef,Nspan,8);
    fill!(difference,NaN);
    difference[:,1:6] = solutione[:,1:6] - NotInfsolutione[:,1:6]
    difference[:,7:8] = controlse[:,1:2] - NotInfcontrolse[:,1:2]

    #States
    fig, estados_e=plt.subplots(4,2)
    fig.suptitle("Optimal States - Entrepreneurs Problem")
        #ue:
    estados_e[1,1].plot(θespan[:], difference[:,1])
    estados_e[1,1].set(ylabel="ue inf - ue without inf", xlabel = "θe")
        #μe:
    estados_e[1,2].plot(θespan[:], difference[:,2])
    estados_e[1,2].set(ylabel="μe inf - μe without inf", xlabel = "θe")
        #Ye:
    estados_e[2,1].plot(θespan[:], difference[:,3])
    estados_e[2,1].set(ylabel="Ye inf - Ye without inf", xlabel = "θe")
        #λe:
    estados_e[2,2].plot(θespan[:], difference[:,4])
    estados_e[2,2].set(ylabel="λe inf - λe without inf", xlabel = "θe")
        #Lfe:
    estados_e[3,1].plot(θespan[:], difference[:,5])
    estados_e[3,1].set(ylabel="Lfe inf - λe without inf", xlabel = "θe")
        #ωf:
    estados_e[3,2].plot(θespan[:], difference[:,6])
    estados_e[3,2].set(ylabel="ωfe inf - ωfe without inf", xlabel = "θe")
        #z_e:
    estados_e[4,1].plot(θespan[:], difference[:,7])
    estados_e[4,1].set(ylabel="ze inf - ze without inf", xlabel = "θe")
        #n_e:
    estados_e[4,2].plot(θespan[:], difference[:,8])
    estados_e[4,2].set(ylabel="ne inf - ne without inf", xlabel = "θe")







#Global Problem (Reverse):
    uw_end     = solutione[1,1]
    μ_end      = solutione[1,2]
    e_end      = elb
    he_ub      = pa.he(pa.θ_w_ub,e_end)
    ϕ_e_end    = -(( pa.indicator*solutione[1,1]^pa.ϕ +solutione[1,4]*(e_end*controlse[1,2]^pa.α- pa.δ/(1.0+pa.γ)*controlse[1,3]^(1.0+pa.γ)
                 - pa.β/(1.0+pa.σ)*controlse[1,1]^(1.0+pa.σ) - solutione[1,1]) - solutione[1,6]*(controlse[1,2]-controlse[1,3])
                 - solutione[1,8]*controlse[1,3])*he_ub  + solutione[1,2]*controlse[1,2]^pa.α*(1.0 - pa.β*controlse[1,1]^pa.σ))
    y_agg_end  = solutione[1,3]
    λ_end      = 1.0
    lf_agg_end  = solutione[1,5]
    ωf_end     = solutione[1,6]
    li_agg_end = solutione[1,7]
    ωi_end     = solutione[1,8]
    wi_end     = solutione[1,9]
    ϕw_end     = solutione[1,10]
    lf_new_end  = solutione[1,11]
    li_new_end = solutione[1,12]
    y_new_end  = solutione[1,13]
    lfstar_end  = 0.0
    listar_end = 0.0
    ystar_end  = 0.0



    Nspan = 500
    y_end= [uw_end, μ_end, e_end, ϕ_e_end, y_agg_end, λ_end, lf_agg_end, ω_end, lf_new_end, y_new_end, 0.0, 0.0];
    xlb= pa.θ_w_lb;
    xub = pa.θ_w_ub;
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xub:-xstep:xlb;
    solution = Array{Float64}(undef,Nspan,12);
    fill!(solution,NaN);
    my_runge_kutta_reverse!(solution,y_end,xspan,xstep,pa)

    #using DelimitedFiles
    #writedlm("SolutionNew.csv",solution,';')

    θspan = Array{Float64,1}
    θspan = collect(xlb:xstep:xub)

    controls = Array{Float64}(undef,Nspan,4);
    fill!(controls,NaN);
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
        fill!(τ_prime,NaN);
        marginal_taxes(controls, θspan, solution, pa)

        (taxes, taxes_ent, taxes_rev, taxes_rev_ent) = taxes_path(controls, θspan, solution, pa, θespan, solutione, controlse, τ_prime, τ_prime_e);

        taxes_full     = fill(NaN,Nspan,6);
        taxes_ent_full = fill(NaN,Nspan,4);
        taxes_full[:,1:3] = taxes;
        taxes_full[:,4:6] = taxes_rev;
        taxes_ent_full[:,1:2] = taxes_ent;
        taxes_ent_full[:,3:4] = taxes_rev_ent;

        #taxes = DataFrame(τ_prime)
        #names!(taxes,[:tau_c, :tau_n, :tau_l] )
        #taxes2=hcat(DataFrame(theta=θspan),taxes)
        #CSV.write("marginal_taxes.csv",taxes2)

        #Utilities:
        utilities_prime = Array{Float64}(undef,Nspan,2);
        fill!(utilities_prime,NaN);

        for i=1:Nspan
            utilities_prime[i,1] = pa.χ*controls[i,3]^(1.0+pa.ψ)/θspan[i]; #u_w
            utilities_prime[i,2] = controlse[i,2]^pa.α*(1.0-pa.β*controlse[i,1]^pa.σ); #u_e
        end

        #A:
        A_matrix = Array{Float64}(undef,Nspan,3);
        fill!(A_matrix,NaN);

        for i=1:Nspan
            A_matrix[i,1] = (pa.indicator*solution[i,1]^pa.ϕ-solution[i,6]*solution[i,1]) + solution[i,4]/pa.he(θspan[i],solution[i,3]); #A
            n_full_info   = ((solution[i,6]*pa.α*solution[i,3])/solution[i,8])^(1.0/(1.0-pa.α));
            A_matrix[i,2] = -(1.0-pa.α)/pa.α*solution[i,8]*n_full_info; #A if z is 0 and n_full_info
            A_matrix[i,3] = controls[i,2]^pa.α*(1.0-τ_prime[i,1]); #u_e in global
        end

        #Max evasion:
        mat_for_z = Array{Float64}(undef,Nspan,4);
        fill!(mat_for_z,NaN);

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

#ctrlvec = controls
#θvec = θspan
#solvec = solution
#controls,θspan,solution,τ_prime,pa,θespan,solutione

        #Plots:
        graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub,
                bound_e,τ_prime,τ_prime_e,"C:\\Users\\mariagon\\Documents\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs\\Graphs",
                utilities_prime,A_matrix,mat_for_z, proposition1, proposition2, proposition3, taxes_full, taxes_ent_full)
        graphs!(solution,solutione,controls,controlse, θspan, θespan, pa.θ_w_ub,
                bound_e,τ_prime,τ_prime_e,"C:\\Users\\marya\\Documents\\GitHub\\OptimalTaxation\\Planner\\GlobalAndEntrepreneurs\\Graphs",
                utilities_prime,A_matrix,mat_for_z, proposition1, proposition2, proposition3, taxes_rev, taxes_rev_ent)



                diff = fill(NaN,Nspan,2)
                for j = 1:Nspan
                    diff[j,1] = proposition3[j,1]+ proposition3[j,2]-proposition3[j,3]
                    #diff[j,2] = proposition1[j,1]- proposition1[j,2]-proposition1[j,3]
                    diff[j,2] = -solution[j,2]- proposition1[j,2]-proposition1[j,3]
                end

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
                                    prop1[1].plot(θspan[:], -solution[:,2])
                                    prop1[1].plot(θspan[:], proposition1[:,2])
                                    prop1[1].plot(θspan[:], proposition1[:,3])
                                    prop1[1].legend(["-μ","Tl'/(1-Tl') εl/(1+εl) θw hw","εz z/n^α he"],loc="upper right")
                                    prop1[2].plot(θspan[:], diff[:,2])

                                    diff()

        to_plot = fill(NaN,Nspan,3);
        for j=1:Nspan
            #Definition of states and controls:
            θ = θspan[j];

                  uw    = solution[j,1];
                  μ     = solution[j,2];
                  e     = solution[j,3];
                  ϕ_e   = solution[j,4];
                  y_agg = solution[j,5];
                  λ     = solution[j,6];
                  l_agg = solution[j,7];
                  ω     = solution[j,8];

                  zz    = controls[j,1];
                  nn    = controls[j,2];
                  ll    = controls[j,3];
                  pp    = controls[j,4];

                  Tn = taxes_full[j,2];

                  #Find the marginal taxes:
                  to_plot[j,1] = e*nn^pa.α-ω*nn
                  to_plot[j,2] = e*nn^pa.α-ω*nn - Tn
                  to_plot[j,3] = e*nn^pa.α-ω*nn - Tn - zz

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
                    tax_ent[2,1].plot(solution[:,3], to_plot[:,1])
                    tax_ent[2,1].set(ylabel="e*nn^pa.α-ω*e*nn", xlabel="e")
                        #τ_n:
                    tax_ent[2,2].plot(solution[:,3], to_plot[:,2])
                    tax_ent[2,2].set(ylabel="e*nn^pa.α-ω*e*nn - Tn", xlabel="e")
                        #τ_c:
                    tax_ent[2,3].plot(solution[:,3], to_plot[:,3])
                    tax_ent[2,3].set(ylabel="e*nn^pa.α-ω*e*nn - Tn - z", xlabel="e")



                    graphs = fill(NaN,Nspan,2)
                    for j = 1:Nspan
                        graphs[j,1] = contols[:,2]^pa.α


                    end
                    graphs[2:500,2] = diff(controls[:,1])

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
                                        tax_ent[2,1].plot(solution[:,3], to_plot[:,1])
                                        tax_ent[2,1].set(ylabel="e*nn^pa.α-ω*e*nn", xlabel="e")
                                            #τ_n:
                                        tax_ent[2,2].plot(solution[:,3], to_plot[:,2])
                                        tax_ent[2,2].set(ylabel="e*nn^pa.α-ω*e*nn - Tn", xlabel="e")
                                            #τ_c:
                                        tax_ent[2,3].plot(solution[:,3], to_plot[:,3])
                                        tax_ent[2,3].set(ylabel="e*nn^pa.α-ω*e*nn - Tn - z", xlabel="e")
