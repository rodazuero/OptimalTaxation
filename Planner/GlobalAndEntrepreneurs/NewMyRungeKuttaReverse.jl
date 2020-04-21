
function my_runge_kutta_reverse!(solution::Array{Float64},y_end,xspan,step,pa,alg; verbose = false)

    println("Solving differencial equations with RK package from Julia.")
    
    #θw is the upper bound of workers distribution:
    solution[end,:] = y_end;
    (Nspan,) =  size(xspan);

    #The initial value for the auxiliary states:
    solution[Nspan,11] = solution[Nspan,7]./solution[Nspan,9];
    solution[Nspan,12] = solution[Nspan,5]./solution[Nspan,10];

    #Loop over values of θw:
    for i = Nspan:-1:2

        println("i = ", i) #Actual θw

        #Current value for θw:
        x   = xspan[i];
        θ   = xspan[i];
        θ_1 = xspan[i-1]; #Previous θw (we are iterating backwards).

        #The initial values of the states:
        uw0    = solution[i,1];
        μ0     = solution[i,2];
        e0     = solution[i,3];
        ϕ_e0   = solution[i,4];
        y_agg0 = solution[i,5];
        λ0     = solution[i,6];
        l_agg0 = solution[i,7];
        ω0     = solution[i,8];
        l_new0 = solution[i,9];
        y_new0 = solution[i,10];

        initial_val = [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0, l_new0, y_new0];
        θ_limits    = (θ,θ_1); #The limits to solve the differencial equation.

        prob  = ODEProblem(find_states!,initial_val,θ_limits,pa)
        #solution_RKPack = solve(prob)
        solution_RKPack = solve(prob,alg)
        #solution_RKPack = solve(prob,alg_hints=[:stiff])

        #Update vector of states:
        solution[i-1,1:10] = solution_RKPack(θ_1);
        solution[i-1,11]   = solution[i-1,7]./solution[i-1,9];
        solution[i-1,12]   = solution[i-1,5]./solution[i-1,10];

    end
end
