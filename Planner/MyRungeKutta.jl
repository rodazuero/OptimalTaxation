
function my_runge_kutta!(solution::Array{Float64},y0,xspan,step,pa; verbose = false)

    #
    solution[1,:] = y0
    #Allocate memory. I follow the notation in Juddm page 345.
    z1 = Array{Float64,1}(undef,10)
    z2 = Array{Float64,1}(undef,10)
    z3 = Array{Float64,1}(undef,10)
    z4 = Array{Float64,1}(undef,10)

    y = Array{Float64,1}(undef,10)
    ini = Array{Float64,1}(undef,11)

    (Nspan,) =  size(xspan)


    #Loop over values of theta

   for i=1:Nspan
        #Current value for theta
        x = xspan[i];
        #θ = exp(xspan[i]);
        θ = xspan[i]
        #println("it = ", i, " x = ", x, " theta = ", θ)

        #Convert state object into a vector
        y[1] = solution[i,1]  # uw
        y[2] = solution[i,2]  # μ
        y[3] = solution[i,3]  # e;
        y[4] = solution[i,4]  # ϕ_e;
        y[5] = solution[i,5]  # y_agg;
        y[6] = solution[i,6]  # λ;
        y[7] = solution[i,7]  # l_agg;
        y[8] = solution[i,8]  # ω;
        y[9] = solution[i,9]  # l_new;
        y[10] = solution[i,10]  # y_new;

        ini[1] = solution[i,1]  # uw
        ini[2] = solution[i,2]  # μ
        ini[3] = solution[i,3]  # e;
        ini[4] = solution[i,4]  # ϕ_e;
        ini[5] = solution[i,5]  # y_agg;
        ini[6] = solution[i,6]  # λ;
        ini[7] = solution[i,7]  # l_agg;
        ini[8] = solution[i,8]  # ω;
        ini[9] = solution[i,9]  # l_new;
        ini[10] = solution[i,10]  # y_new;
        ini[11] = θ  #Actual \theta_w;

        position=i
        #println(" uw = ", y[1], " mu = ", y[2], " e = ", y[3], " phie = ", y[4])
        #println(" Y = ", y[5], " lambda = ", y[6], " L = ", y[7], " omega = ", y[8])
        #We save the approximation of the derivative of the states, according to the order 4 Runge-Kutta method
        find_states!(z1, y, pa, θ, ini);
        any(isnan,z1) && error("z1 is NaN ")
        find_states!(z2, y+0.5*step*z1, pa, x+0.5*step, ini );
        any(isnan,z2) && error("z2 is NaN ")
        find_states!(z3, y+0.5*step*z2, pa, x+0.5*step, ini );
        any(isnan,z3) && error("z3 is NaN ")
        find_states!(z4, y+step*z3, pa, x+step, ini);
        any(isnan,z4) && error("z4 is NaN ")

        if i==1
            #=agg[i,1]=z1[7]
            agg[i,2]=z1[5]
            agg[i,3]=z1[9]
            agg[i,4]=z1[10]
            solution[i,11]=agg[i,1]./agg[i,3] #L*
            solution[i,12]=agg[i,3]./agg[i,4] #Y*=#

            solution[i,11]=z1[7]./z1[9]
            solution[i,12]=z1[5]./z1[10]

        end


        #Update vector of states
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;
        if i < Nspan
            solution[i+1,1:10] = solution[i,1:10] + step*dy;
            #println(i,θ,"   ", solution[i,:])
            #=agg[i+1,1]=solution[i+1,7]
            agg[i+1,2]=solution[i+1,5]
            agg[i+1,3]=solution[i+1,9]
            agg[i+1,4]=solution[i+1,10]
            solution[i+1,11]=agg[i+1,1]./agg[i+1,3]
            solution[i+1,12]=agg[i+1,3]./agg[i+1,4]=#
            solution[i+1,11]=solution[i+1,7]./solution[i+1,9]
            solution[i+1,12]=solution[i+1,5]./solution[i+1,10]
        end


        #solution[i,11:12]=[0,0]

        verbose && println(" duw = ", dy[1], " dmu = ", dy[2], " de = ", dy[3], " dphie = ", dy[4])
        verbose && println(" dY = ", dy[5], " dlambda = ", dy[6], " dL = ", dy[7], " domega = ", dy[8])
        verbose && println("    ")

    end
end
