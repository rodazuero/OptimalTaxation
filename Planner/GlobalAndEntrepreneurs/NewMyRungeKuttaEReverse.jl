function my_runge_kuttae_reverse!(statese::Array{Float64,2},xspan::Array{Float64,1},y_end::Array{Float64,1}, step, pa, alg; verbose::Bool = true)

    verbose && println("Solving differencial equations with RK package from Julia.")

    # θe is the upper bound of entrepreneurs distribution:
    statese[end,1:8] = y_end;
    (Nspan,) =  size(xspan);

    #The initial value for the auxiliary states:
    statese[end,9]  = statese[end,5]/statese[end,8];
    statese[end,10] = statese[end,3]/statese[end,7];

    for i=Nspan:-1:2

        verbose && println("i = ", i) #Actual θe

        #Current value for θe
        θe   = xspan[i]
        θe_1 = xspan[i-1] #Previous θe (we are iterating backwards)

        #The values of the initial point:
        ue0     = statese[i,1];
        μe0     = statese[i,2];
        ye_agg0 = statese[i,3];
        λe0     = statese[i,4];
        le_agg0 = statese[i,5];
        ωe0     = statese[i,6];
        le_new0 = statese[i,7];
        ye_new0 = statese[i,8];

        initial_val = [ue0, μe0, ye_agg0, λe0, le_agg0, ωe0,le_new0,ye_new0];
        θe_limits   = (θe,θe_1); #The limits to solve the differencial equation.
        #println("θe_limits = ", θe_limits)

        prob = ODEProblem(find_statese!,initial_val,θe_limits,pa);
        #solution_RKPack = solve(prob)
        solution_RKPack = solve(prob,alg, save_everystep=false)
        #solution_RKPack = solve(prob,alg_hints=[:stiff])

        #Update vector of states:
        statese[i-1,1:8] = solution_RKPack(θe_1);
        statese[i-1,9]   = statese[i-1,5]./statese[i-1,7];
        statese[i-1,10]  = statese[i-1,3]./statese[i-1,8];

    end
end
