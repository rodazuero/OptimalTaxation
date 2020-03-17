function my_runge_kuttae_reverse!(solution::Array{Float64},y_end,xspan,step,pa,alg; verbose = false)

    # θe is the upper bound of entrepreneurs distribution:
    solution[end,1:8] = y_end;
    (Nspan,) =  size(xspan);

    #The initial value for the auxiliary states:
    solution[end,9]  = solution[end,5]/solution[end,8];
    solution[end,10] = solution[end,3]/solution[end,7];

    #Loop over values of θ_e:
    for i=Nspan:-1:2

        #println("i = ", i) #Actual θe

        #Current value for θ_e
        θe   = xspan[i]
        θe_1 = xspan[i-1]

        #The values of the initial point:
        ue0     = solution[i,1];
        μe0     = solution[i,2];
        ye_agg0 = solution[i,3];
        λe0     = solution[i,4];
        le_agg0 = solution[i,5];
        ωe0     = solution[i,6];
        le_new0 = solution[i,7];
        ye_new0 = solution[i,8];

        initial_val = [ue0, μe0, ye_agg0, λe0, le_agg0, ωe0,le_new0,ye_new0];
        θe_limits   = (θe,θe_1); #The limits to solve the differencial equation.

        prob = ODEProblem(find_statese!,initial_val,θe_limits,pa);
        #solution_RKPack = solve(prob)
        solution_RKPack = solve(prob,alg)
        #solution_RKPack = solve(prob,alg_hints=[:stiff])

        #Update vector of states:
        solution[i-1,1:8] = solution_RKPack(θe_1);
        solution[i-1,9]   = solution[i-1,5]./solution[i-1,7];
        solution[i-1,10]  = solution[i-1,3]./solution[i-1,8];

    end
end
