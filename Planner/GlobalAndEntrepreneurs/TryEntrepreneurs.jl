
function my_runge_kuttae_reverse!(solution::Array{Float64},y_end,xspan,step,pa,θw; verbose = false)

    # θw is the upper bound of workers distribution. It´s taken from the global problem.
    solution[end,1:8] = y_end;
    solution[end,9]  = solution[end,5]/solution[end,8];
    solution[end,10] = solution[end,3]/solution[end,7];

    #Allocate memory. I follow the notation in Juddm page 345.
    z1 = Array{Float64,1}(undef,8);
    z2 = Array{Float64,1}(undef,8);
    z3 = Array{Float64,1}(undef,8);
    z4 = Array{Float64,1}(undef,8);

    y   = Array{Float64,1}(undef,8);
    ini = Array{Float64,1}(undef,9);

    (Nspan,) =  size(xspan);

    #Loop over values of θ_e:
    for i=Nspan:-1:1

        println("i = ", i)
        #Current value for θ_e
        x = xspan[i];
        #θ = exp(xspan[i]);
        θe = xspan[i]
        #println("it = ", i, " x = ", x, " theta_e = ", θe)

        #Convert state object into a vector
        for j = 1:8
            y[j]   = solution[i,j]
            ini[j] = solution[i,j]
        end
            ini[9] = θe  #Actual θ_e;
            #println("θeRK = ", ini[end])

        #println(" ue = ", y[1], " μe = ", y[2], " ye = ", y[3], " λe = ", y[4])
        #println(" Le = ", y[5], " ωe = ", y[6], " Le_new = ", y[7], " Ye_new = ", y[8])

        #We save the approximation of the derivative of the states, according to the order 4 Runge-Kutta method
        #println("z1")
        find_statese!(z1, y, pa, θw, θe, ini);
        any(isnan,z1) && error("z1 is NaN ")
        #println("z2")
        find_statese!(z2, y+0.5*step*z1, pa, θw, x+0.5*step, ini );
        any(isnan,z2) && error("z2 is NaN ")
        #println("z3")
        find_statese!(z3, y+0.5*step*z2, pa, θw, x+0.5*step, ini );
        any(isnan,z3) && error("z3 is NaN ")
        #println("z4")
        find_statese!(z4, y+step*z3, pa, θw, x+step, ini);
        any(isnan,z4) && error("z4 is NaN ")

        #Update vector of states
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;
        if i > 1
            solution[i-1,1:8] = solution[i,1:8] - step*dy;
            solution[i-1,9]   = solution[i-1,5]./solution[i-1,7]
            solution[i-1,10]  = solution[i-1,3]./solution[i-1,8]
        end
    end
end
