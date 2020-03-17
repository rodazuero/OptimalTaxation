function my_runge_kutta_reverse!(solution::Array{Float64},y_end,xspan,step,pa; verbose = false)
    #
    solution[end,:] = y_end
    #agg[1,:] = [0,0,0,0]
    #Allocate memory. I follow the notation in Juddm page 345.
    z1 = Array{Float64,1}(undef,10)
    z2 = Array{Float64,1}(undef,10)
    z3 = Array{Float64,1}(undef,10)
    z4 = Array{Float64,1}(undef,10)

    y = Array{Float64,1}(undef,10)
    ini = Array{Float64,1}(undef,11)

    (Nspan,) =  size(xspan)


    #Loop over values of θw

    for i = 1:Nspan

        #println("i = ", i)
        #Current value for theta
        x = xspan[i];
        #θ = exp(xspan[i]);
        θ = xspan[i]
        #println("it = ", i, " x = ", x, " theta = ", θ)

        #Convert state object into a vector
        y[1] = solution[Nspan+1-i,1]  # uw
        y[2] = solution[Nspan+1-i,2]  # μ
        y[3] = solution[Nspan+1-i,3]  # e;
        y[4] = solution[Nspan+1-i,4]  # ϕ_e;
        y[5] = solution[Nspan+1-i,5]  # y_agg;
        y[6] = solution[Nspan+1-i,6]  # λ;
        y[7] = solution[Nspan+1-i,7]  # l_agg;
        y[8] = solution[Nspan+1-i,8]  # ω;
        y[9] = solution[Nspan+1-i,9]  # l_new;
        y[10] = solution[Nspan+1-i,10]  # y_new;

        ini[1] = solution[Nspan+1-i,1]  # uw
        ini[2] = solution[Nspan+1-i,2]  # μ
        ini[3] = solution[Nspan+1-i,3]  # e;
        ini[4] = solution[Nspan+1-i,4]  # ϕ_e;
        ini[5] = solution[Nspan+1-i,5]  # y_agg;
        ini[6] = solution[Nspan+1-i,6]  # λ;
        ini[7] = solution[Nspan+1-i,7]  # l_agg;
        ini[8] = solution[Nspan+1-i,8]  # ω;
        ini[9] = solution[Nspan+1-i,9]  # l_new;
        ini[10] = solution[Nspan+1-i,10]  # y_new;
        ini[11] = θ  #Actual \theta_w;

        #println(" uw = ", y[1], " mu = ", y[2], " e = ", y[3], " phie = ", y[4])
        #println(" Y = ", y[5], " lambda = ", y[6], " L = ", y[7], " omega = ", y[8])
        #We save the approximation of the derivative of the states, according to the order 4 Runge-Kutta method
        #println("z1")

        A_cons = (pa.indicator*y[1]^pa.ϕ-y[6]*y[1])+y[4]/pa.he(ini[11],y[3]);
        #println("States_A_cons = ", A_cons)

        find_states!(z1, y, pa, θ, ini);
        any(isnan,z1) && error("z1 is NaN ")

        #println("z2")
        find_states!(z2, y+0.5*step*z1, pa, x+0.5*step, ini );
        any(isnan,z2) && error("z2 is NaN ")

        #println("z3")
        find_states!(z3, y+0.5*step*z2, pa, x+0.5*step, ini );
        any(isnan,z3) && error("z3 is NaN ")

        #println("z4")
        find_states!(z4, y+step*z3, pa, x+step, ini);
        any(isnan,z4) && error("z4 is NaN ")

        if i==1
            solution[Nspan,11]= solution[Nspan,7]./solution[Nspan,9]
            solution[Nspan,12]=solution[Nspan,5]./solution[Nspan,10]
        end

        #Update vector of states
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;

        if i < Nspan
            solution[Nspan-i,1:10] = solution[Nspan+1-i,1:10] - step*dy;
            solution[Nspan-i,11]=solution[Nspan-i,7]./solution[Nspan-i,9]
            solution[Nspan-i,12]=solution[Nspan-i,5]./solution[Nspan-i,10]
        end
    end

end
