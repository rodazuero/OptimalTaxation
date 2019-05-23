function find_uw0(y0, uw0_low, uw0_up, Nspan, pa)

    uw0guess    = y0[1];
    μ0     = y0[2];
    e0     = y0[3];
    ϕ_e0   = y0[4];
    y_agg0 = y0[5];
    λ0     = y0[6];
    l_agg0 = y0[7];
    ω0     = y0[8];

    final_mu_aux(uw)= final_mu([uw, μ0, e0, ϕ_e0, y_agg0, λ0, l_agg0, ω0] ,Nspan, pa)
    u0 = find_zero(final_mu_aux, (uw0_low, uw0_up))
end



function final_mu(y0::Array{Float64,1},Nspan::Int64, pa::Param)
    #log-scale for the range of theta,
    xlb= log(pa.θ_w_lb);
    xub = log(pa.θ_w_ub);
    xstep = (xub - xlb)/(Nspan - 1);
    xspan = xlb:xstep:xub;
    solution = Array{Float64}(undef,Nspan,8);
#println("hello")
    try
        #Solve the problem, return value of mu at theta_ub
        my_runge_kutta!(solution,y0,xspan,xstep,pa)
        return μ_ub = solution[end,2]
    catch
        # When no solution, mu is positive, so I return any posutive value, e.g 1
        return μ_ub = 1.0
    end
end
