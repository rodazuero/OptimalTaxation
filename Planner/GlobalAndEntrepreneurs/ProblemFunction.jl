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


function graphs!(solution::Array{Float64},solutione::Array{Float64},controls::Array{Float64},controlse::Array{Float64}, θspan, θespan, thetaw_ub, bound_e,τ_prime::Array{Float64},τ_prime_e::Array{Float64},dir::AbstractString,utilities_prime::Array{Float64},A_mat::Array{Float64})

        rc("font", family="serif")
        original_dir=pwd()

        #Directory to save graphs
        cd(dir)

        #Global problem:

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
            fig, estados_aux=plt.subplots(1,2)
            fig.suptitle("Auxiliar States - Global Problem")
                #L*:
            estados_aux[1,1].plot(θspan[1:500], solution[1:500,11])
            estados_aux[1,1].plot(θspan[1:500], repeat([0],500), "tab:green")
            #estados_aux[1,1].set_title("L*")
            estados_aux[1,1].set(ylabel="L*")
                #Y*:
            estados_aux[2].plot(θspan[1:500], solution[1:500,12])
            estados_aux[2].plot(θspan[1:500], repeat([0.15],500), "tab:green")
            #estados_aux[2].set_title("Y*")
            estados_aux[2].set(ylabel="Y*")

            #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\StatesAux.png")
            #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\StatesAux.png")
            savefig("StatesAux.png")

            #Controls
            fig, controles=plt.subplots(2,2)
            fig.suptitle("Optimal Controls - Global Problem")
                #z:
            controles[1,1].plot(θspan[1:500], controls[:,1])
            #controles[1,1].set_title("z")
            controles[1,1].set(ylabel="z")
                #n:
            controles[1,2].plot(θspan[1:500], controls[1:500,2])
            #controles[1,2].set_title("n")
            controles[1,2].set(ylabel="n")
                #l:
            controles[2,1].plot(θspan[1:500], controls[1:500,3])
            #controles[2,1].set_title("l")
            controles[2,1].set(ylabel="l")
                #p:
            controles[2,2].plot(θspan[1:500], controls[1:500,4])
            #controles[2,2].set_title("p")
            controles[2,2].set(ylabel="p")

            #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\Controls.png")
            #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\Controls.png")
            savefig("Controls.png")

            #Marginal taxes:
            fig, margtax=plt.subplots(1,3)
            fig.suptitle("Marginal Taxes")
                #τ_c_prime:
            margtax[1,1].plot(θspan[1:500], τ_prime[:,1])
            #margtax[1,1].set_title("τ_c'")
            margtax[1,1].set(ylabel="τ_c'")
                #τ_n_prime:
            margtax[2].plot(θspan[1:500], τ_prime[:,2])
            #margtax[2].set_title("τ_n'")
            margtax[2].set(ylabel="τ_n'")
                #τ_l_prime:
            margtax[3].plot(θspan[1:500], τ_prime[:,3])
            #margtax[3].set_title("τ_c'")
            margtax[3].set(ylabel="τ_l'")

            #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\Taxes.png")
            #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\Taxes.png")
            savefig("Taxes.png")

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

        #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\StatesEntrepreneurs.png")
        #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\StatesEntrepreneurs.png")
        savefig("StatesEntrepreneurs.png")

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

        #savefig("C:\\Users\\marya\\Dropbox\\OptimalTaxation\\PlannerMaryan\\Results\\StatesAuxEntrepreneurs.png")
        #savefig("C:\\Users\\mariagon\\Dropbox\\Results\\StatesAuxEntrepreneurs.png")
        savefig("StatesAuxEntrepreneurs.png")

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
        savefig("TaxesEntrepreneurs.png")

        #Graphs for the utilities:
        fig, utilities=plt.subplots(1,2)
        fig.suptitle("Change in Utilities")
            #uw_prime:
        utilities[1].plot(θspan[1:500], utilities_prime[:,1])
        utilities[1].set(ylabel="uw'")
            #ue_prime:
        utilities[2].plot(θespan[1:500], utilities_prime[:,2])
        utilities[2].set(ylabel="ue'")

        savefig("ChangeUtilities.png")

        #Graphs for A:
        fig, A_graphs=plt.subplots(1,2)
        fig.suptitle("A")
            #Value for A:
        #A_graphs = plot(θspan[1:500], A_mat[:,1], A_mat[:,2])
        #suptitle("A")
        A_graphs[1] = plot(θspan[1:500], A_mat[:,1])
        A_graphs[1].set(ylabel="A")
            #Bound for A:
        A_graphs[2].plot(θspan[1:500], A_mat[:,2])
        A_graphs[2].set(ylabel="A Bound")

        savefig("AGraphs.png")

        #Return to original directory
        cd(original_dir)
end
