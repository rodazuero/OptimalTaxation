function my_runge_kutta_reverse!(states::Array{Float64,2}, xspan::Array{Float64,1}, y_end::Array{Float64,1}, step, pa, alg, verbose::Bool = false)

    verbose && println("Solving differential equations with RK package from Julia.")

    #θw is the upper bound of workers distribution:
    states[end,:] = y_end;
    (Nspan,) =  size(xspan);

    #The initial value for the auxiliary states:
    states[Nspan,11] = states[Nspan,7]./states[Nspan,9];
    states[Nspan,12] = states[Nspan,5]./states[Nspan,10];

    # No need for loop over values of θw:

    # Bounds for θw:
    θ   = xspan[end];
    θ_1 = xspan[1]; #Previous θw (we are iterating backwards).

    #The initial values of the states:
    uw0    = y_end[1]
    μ0     = y_end[2]
    e0     = y_end[3]
    ϕ_e0   = y_end[4]
    y_agg0 = y_end[5]
    λ0     = y_end[6]
    l_agg0 = y_end[7]
    ω0     = y_end[8]
    l_new0 = y_end[9]
    y_new0 = y_end[10]

    initial_val = [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0, l_new0, y_new0];
    θ_limits    = (θ,θ_1); #The limits to solve the differencial equation.

    prob  = ODEProblem(find_states!,initial_val,θ_limits,pa)
    condition_Θ_e_l(u,t,integrator)=u[3]-pa.θ_e_lb    
    callback_Θ_e_l=ContinuousCallback(condition_Θ_e_l, affect_Θ_e_l!; save_positions=(false,true) )
    solution_RKPack = solve(prob, alg, callback=callback_Θ_e_l; saveat=xspan)
    #solution_RKPack = solve(prob,alg_hints=[:stiff])
    (nstates, lenght_sol)=size(solution_RKPack)
    
    xspan[1:Nspan-lenght_sol].=NaN
    for i = 2:lenght_sol    
        xspan[Nspan+1-i]=solution_RKPack.t[i] 
        #Update vector of states:
        states[Nspan+1-i,1:10] = solution_RKPack[i]
        states[Nspan+1-i,11]   = states[Nspan+1-i,7]./states[Nspan+1-i,9];
        states[Nspan+1-i,12]   = states[Nspan+1-i,5]./states[Nspan+1-i,10];
    end
    return lenght_sol
end

# function condition_Θ_e_l(u,t,integrator)
#     u[3]-6.33745 #pa.θ_e_lb
# end

function affect_Θ_e_l!(integrator)
    println("Hit")
    terminate!(integrator)
end