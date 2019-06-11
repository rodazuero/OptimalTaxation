function find_uw0(y0, uw0_low, uw0_up, Nspan, pa)
    #Give meaningful names to states
    uw0guess    = y0[1];
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
    mu_low*mu_up > 0.0  && error("Change range of uw0")

    #Solve for uw0 such that mu(theta_ub)=0
    uw0 = find_zero(final_mu_aux, (uw0_low, uw0_up))
    #Convergence happens from below, and it is not defined from above
    #hence I need to adjust uw0 downwards a little bit
    uw0_adjusted = 0.9999999999*uw0;
    abs(final_mu_aux(uw0_adjusted)) > 10.0^-3 && error("Final mu is not zero")
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

function final_e(y0::Array{Float64,1}, ϕ_e0_low, ϕ_e0_up ,Nspan::Int64, pa::Param)

    sol = find_uw0(y0, uw0_low, uw0_up, Nspan, pa);



    #Solve the problem, return value of e at theta_ub
    #this value is consitent with mu(theta_ub)=0
    my_runge_kutta!(solution,y0,xspan,xstep,pa)

    solution[end,3] - pa.θ_e_ub, solution[end,2]
end
