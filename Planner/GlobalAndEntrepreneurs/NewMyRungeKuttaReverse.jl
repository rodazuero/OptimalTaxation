
function my_runge_kutta_reverse!(states::Array{Float64,2}, y_end::Array{Float64,1}, xspan, step, pa, alg, verbose::Bool = false)

    verbose && println("Solving differential equations with RK package from Julia.")
    
    #θw is the upper bound of workers distribution:
    states[end,:] = y_end;
    (Nspan,) =  size(xspan);

    #The initial value for the auxiliary states:
    states[Nspan,11] = states[Nspan,7]./states[Nspan,9];
    states[Nspan,12] = states[Nspan,5]./states[Nspan,10];

    #Loop over values of θw:
    for i = Nspan:-1:2

        verbose && println("i = ", i, ) #Actual θw

        # Current value for θw:
        # x   = xspan[i];
        θ   = xspan[i];
        θ_1 = xspan[i-1]; #Previous θw (we are iterating backwards).

        #The initial values of the states:
        uw0    = states[i,1];
        μ0     = states[i,2];
        e0     = states[i,3];
        ϕ_e0   = states[i,4];
        y_agg0 = states[i,5];
        λ0     = states[i,6];
        l_agg0 = states[i,7];
        ω0     = states[i,8];
        l_new0 = states[i,9];
        y_new0 = states[i,10];

        initial_val = [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0, l_new0, y_new0];
        θ_limits    = (θ,θ_1); #The limits to solve the differencial equation.

        prob  = ODEProblem(find_states!,initial_val,θ_limits,pa)
        #solution_RKPack = solve(prob)
        solution_RKPack = solve(prob, alg, save_everystep=false)
        #solution_RKPack = solve(prob,alg_hints=[:stiff])

        #Update vector of states:
        states[i-1,1:10] = solution_RKPack[2];
        states[i-1,11]   = states[i-1,7]./states[i-1,9];
        states[i-1,12]   = states[i-1,5]./states[i-1,10];

    end
end
