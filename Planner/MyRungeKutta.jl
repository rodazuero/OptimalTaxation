
function my_runge_kutta!(solution::Array{Float64},y0,xspan,step)

    #
    solution[1,:] = y0
    #Allocate memory. I follow the notation in Juddm page 345.
    z1 = Array{Float64,1}(undef,8)
    z2 = Array{Float64,1}(undef,8)
    z3 = Array{Float64,1}(undef,8)
    z4 = Array{Float64,1}(undef,8)

    y = Array{Float64,1}(undef,8)

    (Nspan,) =  size(xspan)


    #Loop over values of theta

    for i=1:Nspan
        #Current value for theta
        x = xspan[i];
        θ = exp(xspan[i]);
        println("it = ", i, " x = ", x, " theta = ", θ)

        #Convert state object into a vector
        y[1] = solution[i,1]  # uw
        y[2] = solution[i,2]  # μ
        y[3] = solution[i,3]  # e;
        y[4] = solution[i,4]  # ϕ_e;
        y[5] = solution[i,5]  # y_agg;
        y[6] = solution[i,6]  # λ;
        y[7] = solution[i,7]  # l_agg;
        y[8] = solution[i,8]  # ω;
        println(" uw = ", y[1], " mu = ", y[2], " e = ", y[3], " phie = ", y[4])
        println(" Y = ", y[5], " lambda = ", y[6], " L = ", y[7], " omega = ", y[8])

        #We save the approximation of the derivative of the states, according to the order 4 Runge-Kutta method
        find_states!(z1, y, pa, θ);
        find_states!(z2, y+0.5*step*z1, pa, exp(x+0.5*step) );
        find_states!(z3, y+0.5*step*z2, pa, exp(x+0.5*step) );
        find_states!(z4, y+step*z3, pa, exp(x+0.5*step) );

        #Update vector of states
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;
        if i < Nspan
            solution[i+1,:] = solution[i,:] + step*dy;
        end
        println(" duw = ", dy[1], " dmu = ", dy[2], " de = ", dy[3], " dphie = ", dy[4])
        println(" dY = ", dy[5], " dlambda = ", dy[6], " dL = ", dy[7], " domega = ", dy[8])

        println("    ")

    end
end
