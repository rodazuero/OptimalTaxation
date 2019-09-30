function distance!(x::Vector,grad::Vector, fixed_initial_states, Nspan, pa, solution::Array{Float64}, io::IOStream)
    if length(grad) > 0
        println("unknown gradient")
    end

    io = open("Results.txt", "a");
    #1 Define initial states
    #1.1 fixed initial states
    μ0   = fixed_initial_states[1];
    e0     = fixed_initial_states[2];
    y_agg0 = fixed_initial_states[3];
    l_agg0 = fixed_initial_states[4];
    y_new0 = fixed_initial_states[5];
    l_new0 = fixed_initial_states[6];
    ystar0 = fixed_initial_states[7];
    lstar0 = fixed_initial_states[8];
    #1.2 initial states that are going to be moving for the shooting
    uw0    = x[1];
    ϕ_e0   = x[2];
    λ0     = x[3];
    ω0     = x[4];
    println(" ")
    println("uw0 = $uw0, ϕ_e0 = $ϕ_e0, λ0 = $λ0, ω0= $ω0, ")
    #println(io," ")
    #print(io, "$uw0, $ϕ_e0, $λ0, $ω0, ")

    #2. Solve differential equation given inital states
    #2.1 log-scale for the range of theta,
    xlb= log(pa.θ_w_lb);
    xub = log(pa.θ_w_ub);
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xlb:xstep:xub;
    #2.2 Define vector for initial guess
    y0= [uw0, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0, l_new0, y_new0, ystar0, lstar0]
    try
        #Solve the problem, return the euclidean distance from solution to the boundary conditions
        my_runge_kutta!(solution,y0,xspan,xstep,pa)
        #av_mu = mean(solution[end,2]);
        rel_error_mu = (solution[end,2] - 0.0);
        #av_e = mean(solution[end,3]);
        rel_error_e = (solution[end,3] - pa.θ_e_ub)/pa.θ_e_ub;
        #av_Y = mean(solution[end,5]);
        rel_error_Y = ( solution[end,12] - pa.G) ;
        #av_L = mean(solution[end,7]);
        rel_error_L = ( solution[end,11] - 0.0 );

        distance =   ( rel_error_mu )^2 +( rel_error_e )^2 +( rel_error_Y )^2 +(rel_error_L )^2

        println("mu = $(solution[end,2]) , delta_e = $(solution[end,3] - pa.θ_e_ub), Y = $(solution[end,5]), L= $(solution[end,7]), dist= $distance ")
        #print(io, "$(solution[end,2]) , $(solution[end,3] - pa.θ_e_ub) , $(solution[end,5]), $(solution[end,7]),  $distance \n")
        #close(io)
        return distance
    catch
        # When no solution, mu is positive, so I return any posutive value, e.g 1
        print(io, " , , , , No solution \n")
        close(io)
        return distance= Inf
    end
end



function find_uw0(y0, uw0_low, uw0_up, Nspan, pa)
    #Give meaningful names to states
    uw0guess    = y0[1]
    μ0     = y0[2];
    e0     = y0[3];
    ϕ_e0   = y0[4];
    y_agg0 = y0[5];
    λ0     = y0[6];
    l_agg0 = y0[7];
    ω0     = y0[8];

    #Define handle function to be minimized
    final_mu_aux(uw)= final_mu([uw, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0] ,Nspan, pa)

    #Make sure the function changes it sign in the inteval
    mu_low = final_mu_aux(uw0_low)
    mu_up = final_mu_aux(uw0_up)
    mu_low*mu_up > 0.0  && error("Change range of uw0. mu(uw_low) = ", mu_low, " and mu(uw_up ) = ", mu_up)

    #Solve for uw0 such that mu(theta_ub)=0
    uw0 = find_zero(final_mu_aux, (uw0_low, uw0_up))
    #Convergence happens from below, and it is not defined from above
    #hence I need to adjust uw0 downwards a little bit
    uw0_adjusted = 0.9999999999*uw0;
    #abs(final_mu_aux(uw0_adjusted)) > 10.0^-3 && error("Final mu is not zero")
    y0[1]=uw0_adjusted;

    #Define grid to solve the differential equation
    xlb= log(pa.θ_w_lb);
    xub = log(pa.θ_w_ub);
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xlb:xstep:xub;
    solution = Array{Float64}(undef,Nspan,8);

    #Solve the problem, return value of e at theta_ub
    my_runge_kutta!(solution,y0,xspan,xstep,pa)

    solution
end



function final_mu(y0::Array{Float64,1},Nspan::Int64, pa::Param)
  # This funtion returns mu(theta_ub), given a vector of initial states y0
  # and a number of values of theta to build the grid for Runge-Kutta

    #log-scale for the range of theta,
    xlb= log(pa.θ_w_lb);
    xub = log(pa.θ_w_ub);
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xlb:xstep:xub;
    solution = Array{Float64}(undef,Nspan,8);
    try
        #Solve the problem, return value of mu at theta_ub
        my_runge_kutta!(solution,y0,xspan,xstep,pa)
        return μ_ub = solution[end,2]
    catch
        # When no solution, mu is positive, so I return any posutive value, e.g 1
        return μ_ub = 1.0
    end
end

function final_delta_e(y0::Array{Float64,1}, uw0_low, uw0_up, Nspan::Int64, pa::Param)
    # Given an initial state vector, returns e(theta_w_ub) - theta_e_ub that is consistent with mu(theta_ub)=0
    # Notice that uw0=y0[1] is redundant because it will be replaced by the optimal one.
    solution = find_uw0(y0, uw0_low, uw0_up, Nspan, pa);

    solution[end,3] - pa.θ_e_ub
end



function find_ϕ_e(y0::Array{Float64,1}, ϕ_e0_low, ϕ_e0_up,uw0_low, uw0_up, Nspan::Int64, pa::Param)

    #Give meaningful names to states
    uw0guess  = y0[1];
    μ0     = y0[2];
    e0     = y0[3];
    ϕ_e0_guess = y0[4];
    y_agg0 = y0[5];
    λ0     = y0[6];
    l_agg0 = y0[7];
    ω0     = y0[8];

    #Define handle function to be minimized
    final_e_aux(ϕ_e)= final_delta_e([uw0guess, μ0, e0, ϕ_e, y_agg0, λ0, l_agg0, ω0], uw0_low, uw0_up, Nspan, pa)

    #Make sure the function changes it sign in the inteval
    e_low = final_e_aux(ϕ_e0_low)
    e_up = final_e_aux(ϕ_e0_up)
    e_low*e_up > 0.0  && error("Change range of phie0")

    #Solve for phie0 such that e(theta_w_ub)=theta_e_ub (and uw0 such that mu(theta_ub)=0)
    @time ϕ_e0 = find_zero(final_e_aux, (ϕ_e0_low, ϕ_e0_up));

    #Recover initial value for uw0 (this is very inefficient, but it works)
    y0[4]=1.0000000001*ϕ_e0;
    sol = find_uw0(y0, uw0_low, uw0_up, Nspan, pa);

    #Solve the problem, return value of e at theta_ub
    my_runge_kutta!(solution,y0,xspan,xstep,pa)

end
