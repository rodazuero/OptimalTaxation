function my_runge_kuttae_reverse!(solution::Array{Float64,2},y_end::Array{Float64,1},xspan,step::Float64,pa,θw::Float64,controlsRK::Array{Float64,2}; verbose = false)
#The states vector is given in the following order:
# ue, μe, ye, λe, lfe, ωfe, lie, ωie, wie and ϕwe.
# We get the following auxiliary states: le_new, lie_new and ye_new.

    # θw is the upper bound of workers distribution. It´s taken from the global problem.
    (Nspan,columns) = size(solution);
    num_states      = columns - 3;
    solution[end,1:num_states] = y_end;

    #Defining the vectors for Ruge-Kutta:
    z1  = Array{Float64,1}(undef,num_states);
    z2  = Array{Float64,1}(undef,num_states);
    z3  = Array{Float64,1}(undef,num_states);
    z4  = Array{Float64,1}(undef,num_states);
    y   = Array{Float64,1}(undef,num_states);
    ini = Array{Float64,1}(undef,num_states+1);

    for i ∈ [z1,z2,z3,z4,y,ini]
        fill!(i,NaN);
    end

    #Loop over values of θ_e:
    for i=Nspan:-1:1

        println("i = ", i)
        #Current value for θ_e
        x  = xspan[i];
        θe = xspan[i];
        println(θe)
        #println("it = ", i, " x = ", x, " theta_e = ", θe)

        for j = 1:num_states
            y[j]   = solution[i,j]; #Actual states;
            ini[j] = solution[i,j]; #Actual states;
        end
            ini[end] = θe;  #Actual θe;
            println("θeRK = ", ini[end])

        #We save the approximation of the derivative of the states, according to the order 4 Runge-Kutta method:
        #println("z1")
        find_statese!(z1, y, pa, θw, θe, ini);
        any(isnan,z1) && error("z1 is NaN ")
        #println("z2")
        find_statese!(z2, y+0.5*step*z1, pa, θw, x+0.5*step, ini);
        any(isnan,z2) && error("z2 is NaN ")
        #println("z3")
        find_statese!(z3, y+0.5*step*z2, pa, θw, x+0.5*step, ini);
        any(isnan,z3) && error("z3 is NaN ")
        #println("z4")
        find_statese!(z4, y+step*z3, pa, θw, x+step, ini);
        any(isnan,z4) && error("z4 is NaN ")

        #Update vector of states
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;
        if i > 1
            solution[i-1,1:num_states] = solution[i,1:num_states] - step*dy;

            solution[i-1,14] = solution[i-1,5]./solution[i-1,11];
            solution[i-1,15] = solution[i-1,7]./solution[i-1,12];
            solution[i-1,16] = solution[i-1,3]./solution[i-1,13];
        end

        if i==Nspan
            solution[Nspan,14]=solution[Nspan,5]./solution[Nspan,11];
            solution[Nspan,15]=solution[Nspan,7]./solution[Nspan,12];
            solution[Nspan,16]=solution[Nspan,3]./solution[Nspan,13];
        end
    end
end
