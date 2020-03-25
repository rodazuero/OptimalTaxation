
function my_runge_kuttae_reverse!(solution::Array{Float64},y_end,xspan,step,pa,θw; verbose = false)

    # θw is the upper bound of workers distribution. It´s taken from the global problem.
    solution[end,1:8] = y_end;
    #agg[1,:] = [0,0,0,0]
    #Allocate memory. I follow the notation in Juddm page 345.
    z1 = Array{Float64,1}(undef,8);
    z2 = Array{Float64,1}(undef,8);
    z3 = Array{Float64,1}(undef,8);
    z4 = Array{Float64,1}(undef,8);

    y   = Array{Float64,1}(undef,8);
    ini = Array{Float64,1}(undef,9);

    (Nspan,) =  size(xspan);

    #Loop over values of θ_e:
    for i=1:Nspan

        #println("i = ", i)
        #Current value for θ_e
        x = xspan[i];
        #θ = exp(xspan[i]);
        θe = xspan[i]
        #println("it = ", i, " x = ", x, " theta_e = ", θe)

        #Convert state object into a vector
        y[1] = solution[Nspan+1-i,1]  # ue
        y[2] = solution[Nspan+1-i,2]  # μe
        y[3] = solution[Nspan+1-i,3]  # ye;
        y[4] = solution[Nspan+1-i,4]  # λe;
        y[5] = solution[Nspan+1-i,5]  # le;
        y[6] = solution[Nspan+1-i,6]  # ωe;
        y[7] = solution[Nspan+1-i,7]  # le_new;
        y[8] = solution[Nspan+1-i,8]  # ye_new;

        ini[1] = solution[Nspan+1-i,1]  # ue
        ini[2] = solution[Nspan+1-i,2]  # μe
        ini[3] = solution[Nspan+1-i,3]  # ye;
        ini[4] = solution[Nspan+1-i,4]  # λe;
        ini[5] = solution[Nspan+1-i,5]  # le;
        ini[6] = solution[Nspan+1-i,6]  # ωe;
        ini[7] = solution[Nspan+1-i,7]  # le_new;
        ini[8] = solution[Nspan+1-i,8]  # ye_new;
        ini[9] = θe  #Actual \theta_e;
        println("θeRK = ", ini[end])

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

        if i==1
            #=agg[i,1]=z1[7]
            agg[i,2]=z1[5]
            agg[i,3]=z1[9]
            agg[i,4]=z1[10]
            solution[i,11]=agg[i,1]./agg[i,3] #L*
            solution[i,12]=agg[i,3]./agg[i,4] #Y*=#

            solution[Nspan,9]=solution[Nspan,5]./solution[Nspan,8]
            solution[Nspan,10]=solution[Nspan,3]./solution[Nspan,7]

        end

        #Update vector of states
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;
        if i < Nspan
            solution[Nspan-i,1:8] = solution[Nspan+1-i,1:8] - step*dy;
            solution[Nspan-i,9]=solution[Nspan-i,5]./solution[Nspan-i,7]
            solution[Nspan-i,10]=solution[Nspan-i,3]./solution[Nspan-i,8]
        end

    end


end
