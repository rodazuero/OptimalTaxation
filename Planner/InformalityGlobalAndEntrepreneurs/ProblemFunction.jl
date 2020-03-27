function graphs!(solutione::Array{Float64},NotInfsolutione::Array{Float64},controlse::Array{Float64},NotInfcontrolse::Array{Float64}, θespan::Array{Float64}, τ_prime_e::Array{Float64},NIτ_prime_e::Array{Float64}, difference,dir::AbstractString)

        rc("font", family="serif")
        original_dir=pwd()

        #Directory to save graphs
        cd(dir)

        #States
        fig, estados_e=plt.subplots(5,2)
        fig.suptitle("Optimal States - Entrepreneurs Problem")
            #ue:
        estados_e[1,1].plot(θespan[:], solutione[:,1])
        estados_e[1,1].set_title("u_e")
        estados_e[1,1].set(ylabel="ue")
            #μe:
        estados_e[1,2].plot(θespan[:], solutione[:,2])
        estados_e[1,2].plot(θespan[:], repeat([0],500), "tab:green")
        estados_e[1,2].set_title("μe")
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
        estados_e[3,1].set_title("Le")
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

        savefig("States_Entrepreneurs.png")

        #Controls
        fig, controles_e=plt.subplots(1,3)
        fig.suptitle("Optimal Controls - Entrepreneurs Problem")
            #ze:
        controles_e[1].plot(θespan[:], controlse[:,1])
        controles_e[1].set(ylabel="ze")
            #ne:
        controles_e[2].plot(θespan[:], controlse[:,2])
        controles_e[2].set(ylabel="ne")
            #nie:
        controles_e[3].plot(θespan[:], controlse[:,3])
        controles_e[3].set(ylabel="nie")

        savefig("Controls_Entrepreneurs.png")

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

        savefig("MargTaxes_Entrepreneurs.png")

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

        savefig("CompStates_Entrepreneurs.png")

        #Controls
        fig, controles_e=plt.subplots(1,3)
        fig.suptitle("Optimal Controls - Entrepreneurs Problem")
            #ze:
        controles_e[1].plot(θespan[:], controlse[:,1], θespan[:], NotInfcontrolse[:,1])
        controles_e[1].legend(["ze with informality", "ze without informality"],loc="upper right")
        controles_e[1].set(ylabel="ze", xlabel = "θe")
            #ne:
        controles_e[2].plot(θespan[:], controlse[:,2], θespan[:], NotInfcontrolse[:,2])
        controles_e[2].legend(["ne with informality", "ne without informality"],loc="upper right")
        controles_e[2].set(ylabel="ne", xlabel = "θe")
            #nie:
        controles_e[3].plot(θespan[:], controlse[:,3])
        controles_e[3].set(ylabel="nie", xlabel = "θe")

        savefig("CompControls_Entrepreneurs.png")

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

        savefig("CompMargTaxes_Entrepreneurs.png")

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
        estados_e[3,1].set(ylabel="Lfe inf - Lfe without inf", xlabel = "θe")
            #ωf:
        estados_e[3,2].plot(θespan[:], difference[:,6])
        estados_e[3,2].set(ylabel="ωfe inf - ωfe without inf", xlabel = "θe")
            #z_e:
        estados_e[4,1].plot(θespan[:], difference[:,7])
        estados_e[4,1].set(ylabel="ze inf - ze without inf", xlabel = "θe")
            #n_e:
        estados_e[4,2].plot(θespan[:], difference[:,8])
        estados_e[4,2].set(ylabel="ne inf - ne without inf", xlabel = "θe")

        savefig("CompDifference_Entrepreneurs.png")

        #Return to original directory
        cd(original_dir)
end
