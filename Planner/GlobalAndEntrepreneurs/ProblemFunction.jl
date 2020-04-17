function complete_problem!(solution::Array{Float64},solutione::Array{Float64},controls::Array{Float64},controlse::Array{Float64},y0,xspan,step,pa; verbose = false)

    Nspan = 500

    #Global Problem
    my_runge_kutta!(solution,y0,xspan,xstep,pa)

    θspan = Array{Float64,1}
    θspan = collect(xspan)

    recover_controls!(controls, θspan, solution);

    he_ub= pa.he(pa.θ_w_ub,solution[end,3])
    bound_e= -(pa.indicator*solution[end,1]^pa.ϕ*he_ub+solution[end,6]*(solution[end,3]*controls[end,2]^pa.α -(pa.β/(1.0+pa.σ))*controls[end,1]^(1.0+pa.σ)
            - solution[end,1])*he_ub - solution[end,8]*controls[end,2]*he_ub + solution[end,2]*controls[end,2]^pa.α*(1.0-pa.β*controls[end,1]^pa.σ))
    dist_e= bound_e - solution[end,4]

    #Check adequate end for Y, L and ϕ_e
    #Final L
    solution[end,7] < 0.0 && error("Demand surplus in labor market.", solution[end,7])
    #Final Y
    #solution[end,5]>= 0.0 && error("Supply surplus in goods and services market  ", solution[end,5])
    #Final ϕ_e (ad-hoc value)
    abs(dist_e)> 30.0 && error("Final boundary condition for ϕ_e is not satisfied.", dist_e)

    #Entrepreneurs Problem
    ue0    =   solution[end,1]
    μe0    =   solution[end,2]
    ye0    =   solution[end,5]
    λe0    =   solution[end,6]
    le0    =   solution[end,7]
    ωe0    =   solution[end,8]

    ye0= [ue0, μe0, ye0, λe0, le0, ωe0, 0.0, 0.0];
    elb= solution[end,3];
    eub = pa.θ_e_ub;
    estep = (eub - elb)/(Nspan - 1);
    espan = elb:estep:eub;
    my_runge_kuttae!(solutione,ye0,espan,estep,pa,pa.θ_w_ub)

    θespan = Array{Float64,1}
    θespan = collect(espan)

    recover_controlse!(controlse, pa.θ_w_ub ,θespan, solutione);
    controlse;
end


function graphs!(solution::Array{Float64}, solutione::Array{Float64}, controls::Array{Float64}, controlse::Array{Float64}, θspan::Array{Float64}, θespan::Array{Float64},
                    thetaw_ub, bound_e, τ_prime::Array{Float64}, τ_prime_e::Array{Float64}, dir::AbstractString, utilities_prime::Array{Float64}, A_mat::Array{Float64},
                    mat_z::Array{Float64}, proposition1::Array{Float64}, proposition2::Array{Float64}, proposition3::Array{Float64}, taxes::Array{Float64}, taxes_ent::Array{Float64})

        rc("font", family="serif")
        original_dir=pwd()

        #Directory to save graphs
        cd(dir)

        #Global problem:

            #Full information
            n_full = Array{Float64}(undef,Nspan,1);
            l_full = Array{Float64}(undef,Nspan,1);
            fill!(n_full,NaN);
            fill!(l_full,NaN);
            for i=1:Nspan
                n_full[i]=(pa.α*solution[i,3]/solution[i,8])^(1.0/(1.0-pa.α))
                l_full[i]=((θspan[i]*solution[i,8])/(solution[i,6]*pa.χ))^(1.0/pa.ψ)
            end

            #States
            fig, estados=plt.subplots(4,2)
            fig.suptitle("Optimal States - Global Problem")
                #u_w:
            estados[1,1].plot(θspan[1:500], solution[1:500,1])
            #estados[1,1].set_title("u_w")
            estados[1,1].set(ylabel="u_w")
                #μ:
            estados[1,2].plot(θspan[1:500], solution[1:500,2])
            #estados[1,2].set_title("μ")
            estados[1,2].set(ylabel="μ")
                #e:
            estados[2,1].plot(θspan[1:500], solution[1:500,3])
            estados[2,1].plot(θspan[1:500], repeat([pa.θ_e_ub],500), "tab:green")
            estados[2,1].plot(θspan[1:500], repeat([e_end],500), "tab:blue")
            #estados[2,1].set_title("e")
            estados[2,1].set(ylabel="e")
                #ϕe:
            estados[2,2].plot(θspan[1:500], solution[1:500,4])
            estados[2,2].plot(θspan[1:500], repeat([bound_e],500), "tab:green")
            #estados[2,2].set_title("φe")
            estados[2,2].set(ylabel="φe")
                #Y:
            estados[3,1].plot(θspan[1:500], solution[1:500,5])
            estados[3,1].plot(θspan[1:500], repeat([0.15],500), "tab:green")
            #estados[3,1].set_title("Y")
            estados[3,1].set(ylabel="Y")
                #λ:
            estados[3,2].plot(θspan[1:500], solution[1:500,6])
            #estados[3,2].set_title("λ")
            estados[3,2].set(ylabel="λ")
                #L:
            estados[4,1].plot(θspan[1:500], solution[1:500,7])
            estados[4,1].plot(θspan[1:500], repeat([0],500), "tab:green")
            #estados[4,1].set_title("L")
            estados[4,1].set(ylabel="L")
                #ω:
            estados[4,2].plot(θspan[1:500], solution[1:500,8])
            #estados[4,2].set_title("ω")
            estados[4,2].set(ylabel="ω")

            #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\States.png")
            #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\States.png")
            savefig("States.png")

            #Auxiliar states
            #fig, estados_aux=plt.subplots(1,2)
            #fig.suptitle("Auxiliar States - Global Problem")
                #L*:
            #estados_aux[1,1].plot(θspan[1:500], solution[1:500,11])
            #estados_aux[1,1].plot(θspan[1:500], repeat([0],500), "tab:green")
            #estados_aux[1,1].set_title("L*")
            #estados_aux[1,1].set(ylabel="L*")
                #Y*:
            #estados_aux[2].plot(θspan[1:500], solution[1:500,12])
            #estados_aux[2].plot(θspan[1:500], repeat([0.15],500), "tab:green")
            #estados_aux[2].set_title("Y*")
            #estados_aux[2].set(ylabel="Y*")

            #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\StatesAux.png")
            #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\StatesAux.png")
            #savefig("StatesAux.png")

            #Controls
            fig, controles=plt.subplots(2,2)
            fig.suptitle("Optimal Controls - Global Problem")
                #z:
            controles[1,1].plot(θspan[1:500], controls[:,1])
            #controles[1,1].set_title("z")
            controles[1,1].set(ylabel="z")
                #n:
           controles[1,2].plot(θspan[1:500], controls[1:500,2])
           controles[1,2].plot(θspan[1:500], n_full[1:500], linestyle="--")
           #controles[1,2].set_title("n")
           controles[1,2].set(ylabel="n")
           controles[1,2].legend(["n","n full info"],loc="upper left")
                #l:
            controles[2,1].plot(θspan[1:500], controls[1:500,3])
            controles[2,1].plot(θspan[1:500], l_full[1:500], linestyle="--")
            #controles[2,1].set_title("l")
            controles[2,1].set(ylabel="l")
            controles[2,1].legend(["l","l full info"],loc="upper left")
                #p:
            controles[2,2].plot(θspan[1:500], controls[1:500,4])
            #controles[2,2].set_title("p")
            controles[2,2].set(ylabel="p")

            #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\Controls.png")
            #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\Controls.png")
            savefig("Controls.png")

            #Marginal taxes:
            fig, margtax=plt.subplots(2,3)
            fig.suptitle("Marginal Taxes")
                #τ_c_prime:
            margtax[1,1].plot(θspan[1:500], τ_prime[:,1])
            #margtax[1,1].set_title("τ_c'")
            margtax[1,1].set(ylabel="τ_c'")
                #τ_n_prime:
            margtax[1,2].plot(θspan[1:500], τ_prime[:,2])
            #margtax[2].set_title("τ_n'")
            margtax[1,2].set(ylabel="τ_n'")
                #τ_l_prime:
            margtax[1,3].plot(θspan[1:500], τ_prime[:,3])
            #margtax[3].set_title("τ_c'")
            margtax[1,3].set(ylabel="τ_l'")
                #τ_c:
            margtax[2,1].plot(taxes[:,4], τ_prime[:,1])
            margtax[2,1].set(ylabel="T_c", xlabel="e*n^α- ω n-Tn -z")
                #τ_n:
            margtax[2,2].plot(taxes[:,5], τ_prime[:,2])
            margtax[2,2].set(ylabel="T_n", xlabel="ω n")
                #τ_l:
            margtax[2,3].plot(taxes[:,6], τ_prime[:,3])
            margtax[2,3].set(ylabel="T_l", xlabel="θw l ω")

            #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\Taxes.png")
            #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\Taxes.png")
            savefig("MarginalTaxes.png")

        #Entrepreneurs problem:

        #Full information
           n_full_e = Array{Float64}(undef,Nspan,1);
           for i=1:Nspan
               n_full_e[i]=(pa.α*θespan[i]/solutione[i,6])^(1.0/(1.0-pa.α))
           end

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

        #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\StatesEntrepreneurs.png")
        #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\StatesEntrepreneurs.png")
        savefig("StatesEntrepreneurs.png")

        #Auxiliar states
        # fig, estados_auxE=plt.subplots(1,2)
        # fig.suptitle("Auxiliar States - Global Problem")
        #     #L*:
        # estados_auxE[1,1].plot(θespan[1:500], solutione[1:500,9])
        # estados_auxE[1,1].plot(θespan[1:500], repeat([0],500), "tab:green")
        # #estados_auxE[1,1].set_title("L*")
        # estados_auxE[1,1].set(ylabel="L*")
        #     #Y*:
        # estados_auxE[2].plot(θespan[1:500], solutione[1:500,10])
        # estados_auxE[2].plot(θespan[1:500], repeat([0.15],500), "tab:green")
        # #estados_auxE[2].set_title("Y*")
        # estados_auxE[2].set(ylabel="Y*")

        # #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\StatesAuxEntrepreneurs.png")
        # #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\StatesAuxEntrepreneurs.png")
        # savefig("StatesAuxEntrepreneurs.png")

        #Controls
        fig, controles_e=plt.subplots(1,2)
        fig.suptitle("Optimal Controls - Entrepreneurs Problem")
            #ze:
        controles_e[1,1].plot(θespan[1:500], controlse[:,1])
        #controles_e[1,1].set_title("ze")
        controles_e[1,1].set(ylabel="ze")
            #ne:
        controles_e[2].plot(θespan[1:500], controlse[:,2])
        controles_e[2].plot(θespan[1:500], n_full_e[1:500], linestyle="--")
        #controles_e[1].set_title("ne")
        controles_e[2].set(ylabel="ne")
        controles_e[2].legend(["n","n full info"],loc="upper left")

        #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\ControlsEntrepreneurs.png")
        #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\ControlsEntrepreneurs.png")
        savefig("ControlsEntrepreneurs.png")

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

        #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\TaxesEntrepreneurs.png")
        #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\TaxesEntrepreneurs.png")
        savefig("MarginalTaxesEntrepreneurs.png")

        #Graphs for the utilities:
        fig, utilities=plt.subplots(1,2)
        fig.suptitle("Change in Utilities")
            #ue_prime:
        utilities[1].plot(θspan, utilities_prime[:,1])
        utilities[1].set(ylabel="uw' in global")
        utilities[1].set(xlabel="θ_w")
            #ue_prime:
        utilities[2].plot(solution[:,3], utilities_prime[:,2])
        utilities[2].set(ylabel="ue'")
        utilities[2].set(xlabel="θ_e")

        savefig("MargUtility.png")

        #Graphs for A:
        fig, A_graphs=plt.subplots(1,2)
        fig.suptitle("A Graphs")
            #Value for A:
        #A_graphs = plot(θspan[1:500], A_mat[:,1], A_mat[:,2])
        #suptitle("A")
        A_graphs[1].plot(θspan[1:500], A_mat[:,1])
        A_graphs[1].set(ylabel="A")
            #Bound for A:
        A_graphs[2].plot(θspan[1:500], A_mat[:,2])
        A_graphs[2].set(ylabel="A Bound")

        savefig("AGraphs.png")

        #Graphs for Z:
        fig, Z_graphs=plt.subplots(2,3)
        fig.suptitle("Evasion Graphs")
            #Value for max possible evasion:
            Z_graphs[1,1].plot(θspan[:], mat_z[:,1])
            Z_graphs[1,1].set(ylabel="λen^α - ωn",xlabel="θw")
            #Evasion in model:
            Z_graphs[1,2].plot(θspan[:], mat_z[:,2])
            Z_graphs[1,2].set(ylabel="z",xlabel="θw")
            #Combined:
            Z_graphs[1,3].plot(θspan[:], mat_z[:,1], θspan[:], mat_z[:,2])
            Z_graphs[1,3].set(xlabel="θw")
            Z_graphs[1,3].legend(["λen^α - ωn","z"],loc="upper left")

            #Value for max possible evasion:
            Z_graphs[2,1].plot(θespan[:], mat_z[:,3])
            Z_graphs[2,1].set(ylabel="λθen^α - ωn",xlabel="θe")
            #Evasion in model:
            Z_graphs[2,2].plot(θespan[:], mat_z[:,4])
            Z_graphs[2,2].set(ylabel="z",xlabel="θe")
            #Combined:
            Z_graphs[2,3].plot(θespan[:], mat_z[:,3], θespan[:], mat_z[:,4])
            Z_graphs[2,3].set(ylabel="λen^α - ωn or z", xlabel="θe")
            plt.legend(["λen^α - ωn","z"],loc="upper right")

        savefig("zGraphs.png")

        #Graphs for Propositions:
        #Proposition1:
        #fig, prop1=plt.subplots(1,4)
        #fig.suptitle("Proposition 1")
            #Left Side:
            #prop1[1].plot(θspan[:], -solution[:,2])
            #prop1[1].set(ylabel="-μ",xlabel="θw")
            #Right Side (1):
            #prop1[2].plot(θspan[:], proposition1[:,2])
            #prop1[2].set(ylabel="Tl'/(1-Tl') εl/(1+εl) θw hw", xlabel="θw")
            #Right Side (2):
            #prop1[3].plot(θspan[:], proposition1[:,3])
            #prop1[3].set(ylabel="εz z/n^α he", xlabel="θw")
            #All:
            #prop1[4].plot(θspan[:], -solution[:,2])
            #prop1[4].plot(θspan[:], proposition1[:,2])
            #prop1[4].plot(θspan[:], proposition1[:,3])
            #prop1[4].legend(["-μ","Tl'/(1-Tl') εl/(1+εl) θw hw","εz z/n^α he"],loc="upper right")

            fig, prop1=plt.subplots(1)
            fig.suptitle("Proposition 1")
                prop1.plot(θspan[:], -solution[:,2])
                prop1.plot(θspan[:], proposition1[:,2])
                prop1.plot(θspan[:], proposition1[:,3])
                prop1.legend(["-μ","Tl'/(1-Tl') εl/(1+εl) θw hw","εz z/n^α he"],loc="upper right")
                prop1.set(xlabel="θw")

            savefig("proposition1.png")

        #Proposition2:
        #fig, prop2=plt.subplots(1,3)
        #fig.suptitle("Proposition 2")
            #Left Side:
            #prop2[1].plot(θspan[:], proposition2[:,1])
            #prop2[1].set(ylabel="εz z/n^α",xlabel="θw")
            #Right Side:
            #prop2[2].plot(θspan[:], proposition2[:,2])
            #prop2[2].set(ylabel="e/(1-Tc') Tn'/(1+Tn)", xlabel="θw")
            #All:
            #prop2[3].plot(θspan[:], proposition2[:,1])
            #prop2[3].plot(θspan[:], proposition2[:,2])
            #plt.legend(["εz z/n^α","e/(1-Tc') Tn'/(1+Tn)"],loc="upper right")

            fig, prop2=plt.subplots(1)
            fig.suptitle("Proposition 2")
                prop2.plot(θspan[:], proposition2[:,1])
                prop2.set(xlabel="θw")
                prop2.plot(θspan[:], proposition2[:,2])
                prop2.legend(["εz z/n^α", "e/(1-Tc') Tn'/(1+Tn)"], loc="upper right")

            savefig("proposition2.png")

        #Proposition3:
        #fig, prop3=plt.subplots(1,3)
        #fig.suptitle("Proposition 3")
            #Left Side:
            #prop3[1].plot(θspan[:], proposition3[:,5])
            #prop3[1].set(ylabel="εz z/n^α he + 1/λ [Ve-Vw] g 1/ue'",xlabel="θw")
            #Right Side:
            #prop3[2].plot(θspan[:], proposition3[:,6])
            #prop3[2].set(ylabel="λ he p", xlabel="θw")
            #All:
            #prop3[3].plot(θspan[:], proposition3[:,1])
            #prop3[3].plot(θspan[:], proposition3[:,2])
            #prop3[3].plot(θspan[:], proposition3[:,6])
            #prop3[3].legend(["εz z/n^α he","1/λ [Ve-Vw] g 1/ue'","λ he p"],loc="upper right")

            fig, prop3=plt.subplots(1)
            fig.suptitle("Proposition 3")
                prop3.plot(θspan[:], proposition3[:,1])
                prop3.plot(θspan[:], proposition3[:,2])
                prop3.plot(θspan[:], proposition3[:,3])
                prop3.legend(["εz z/n^α he", "1/λ [Ve-Vw] g 1/ue'", "λ he p"],loc="upper right")
                prop3.set(xlabel="θw")

            savefig("proposition3.png")

            #Graphs for all the span of entrepreneurs (both problems)
            θe_all=vcat(solution[:,3],θespan)

            #States
            fig, estados_e=plt.subplots(3,2)
            fig.suptitle("Optimal States with span of entrepreneurs")
                #u_e:
            estados_e[1,1].plot(θe_all, vcat(solution[1:500,1],solutione[1:500,1]))
            #estados_e[1,1].set_title("u_e")
            estados_e[1,1].set(ylabel="u_e")
                #μ_e:
            estados_e[1,2].plot(θe_all, vcat(solution[1:500,2],solutione[1:500,2]))
            estados_e[1,2].plot(θe_all, repeat([0],1000), "tab:green")
            #estados_e[1,2].set_title("μe")
            estados_e[1,2].set(ylabel="μe")
                #Y_e:
            estados_e[2,1].plot(θe_all, vcat(solution[1:500,5],solutione[1:500,3]))
            estados_e[2,1].plot(θe_all, repeat([0.15],1000), "tab:green")
            #estados_e[2,1].set_title("Ye")
            estados_e[2,1].set(ylabel="Ye")
                #λ_e:
            estados_e[2,2].plot(θe_all, vcat(solution[1:500,6],solutione[1:500,4]))
            #estados_e[2,2].set_title("λe")
            estados_e[2,2].set(ylabel="λe")
                #L_e:
            estados_e[3,1].plot(θe_all, vcat(solution[1:500,7],solutione[1:500,5]))
            estados_e[3,1].plot(θe_all, repeat([0],1000), "tab:green")
            #estados_e[3,1].set_title("Le")
            estados_e[3,1].set(ylabel="Le")
                #ω:
            estados_e[3,2].plot(θe_all, vcat(solution[1:500,8],solutione[1:500,6]))
            #estados_e[3,2].set_title("ωe")
            estados_e[3,2].set(ylabel="ωe")

            savefig("StatesAll.png")

            #Controls
            fig, controles_e=plt.subplots(1,2)
            fig.suptitle("Optimal Controls")
                #ze:
            controles_e[1,1].plot(θe_all, vcat(controls[1:500,1],controlse[1:500,1]))
            #controles_e[1,1].set_title("ze")
            controles_e[1,1].set(ylabel="ze")
                #ne:
            controles_e[2].plot(θe_all, vcat(controls[1:500,2],controlse[1:500,2]))
            controles_e[2].plot(θe_all, vcat(n_full[1:500],n_full_e[1:500]), linestyle="--")
            #controles_e[1].set_title("ne")
            controles_e[2].set(ylabel="ne")
            controles_e[2].legend(["n","n full info"],loc="upper left")

            savefig("ControlsAll.png")

            #Marginal taxes:
            fig, margtax_e=plt.subplots(1,2)
            fig.suptitle("Marginal Taxes")
                #τ_c_prime:
            margtax_e[1].plot(θe_all, vcat(τ_prime[:,1],τ_prime_e[:,1]))
            #margtax_e[1].set_title("τ_c'")
            margtax_e[1].set(ylabel="τ_c'")
                #τ_n_prime:
            margtax_e[2].plot(θe_all, vcat(τ_prime[:,2],τ_prime_e[:,2]))
            #margtax_e[2].set_title("τ_n'")
            margtax_e[2].set(ylabel="τ_n'")

            savefig("MarginalTaxesAll.png")

        #Taxes liabilities:
        fig, tax=plt.subplots(2,3)
        fig.suptitle("Taxes Global Problem")
            #τ_c:
        tax[1,1].plot(solution[:,3], taxes[:,1])
        tax[1,1].set(ylabel="T_c", xlabel="e")
            #τ_n:
        tax[1,2].plot(solution[:,3], taxes[:,2])
        tax[1,2].set(ylabel="T_n", xlabel="e")
            #τ_l:
        tax[1,3].plot(θspan[1:500], taxes[:,3])
        tax[1,3].set(ylabel="T_l", xlabel="θw")
            #τ_c:
        tax[2,1].plot(taxes[:,4], taxes[:,1])
        tax[2,1].set(ylabel="T_c", xlabel="e*n^α- ω n-Tn -z")
            #τ_n:
        tax[2,2].plot(taxes[:,5], taxes[:,2])
        tax[2,2].set(ylabel="T_n", xlabel="ω n")
            #τ_l:
        tax[2,3].plot(taxes[:,6], taxes[:,3])
        tax[2,3].set(ylabel="T_l", xlabel="θw l ω")

        savefig("Taxes_Global.png")

        fig, tax_ent=plt.subplots(2,2)
        fig.suptitle("Taxes Entrepreneurs Problem")
            #τ_c:
        tax_ent[1,1].plot(θespan[1:500], taxes_ent[:,1])
        tax_ent[1,1].set(ylabel="T_c", xlabel="θe")
            #τ_n:
        tax_ent[1,2].plot(θespan[1:500], taxes_ent[:,2])
        tax_ent[1,2].set(ylabel="T_n", xlabel="θe")
            #τ_c:
        tax_ent[2,1].plot(taxes_ent[:,3], taxes_ent[:,1])
        tax_ent[2,1].set(ylabel="T_c", xlabel="θe*n^α- ωe n-Tn-z")
            #τ_n:
        tax_ent[2,2].plot(taxes_ent[:,4], taxes_ent[:,2])
        tax_ent[2,2].set(ylabel="T_n", xlabel="ωe n")

        savefig("Taxes_Entrepreneurs.png")

        #Return to original directory
        cd(original_dir)
end
