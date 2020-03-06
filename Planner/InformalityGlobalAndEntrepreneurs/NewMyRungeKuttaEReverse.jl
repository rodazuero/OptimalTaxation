function my_runge_kuttae_reverse!(solution::Array{Float64,2},y_end::Array{Float64,1},xspan,step::Float64,pa,θw::Float64; verbose = false)
#The states vector is given in the following order:
# ue, μe, ye, λe, le, ωfe, lie, ωie, wie and ϕie.
# We get the following auxiliary states: le_new, lie_new and ye_new.

    # θw is the upper bound of workers distribution. It´s taken from the global problem.
    (Nspan,num_states) = size(solution);
    solution[end,1:num_states] = y_end;

    #Defining the vectors for Ruge-Kutta:
    z1 = Array{Float64,1}(undef,num_states);
    fill!(z1,NaN);
    z2 = Array{Float64,1}(undef,num_states);
    fill!(z2,NaN);
    z3 = Array{Float64,1}(undef,num_states);
    fill!(z3,NaN);
    z4 = Array{Float64,1}(undef,num_states);
    fill!(z4,NaN);

    y   = Array{Float64,1}(undef,num_states);
    fill!(y,NaN);
    ini = Array{Float64,1}(undef,num_states+1);
    fill!(ini,NaN);

    #Loop over values of θ_e:
    for i=Nspan:-1:1

        #println("i = ", i)
        #Current value for θ_e
        x  = xspan[i];
        θe = xspan[i];
        #println("it = ", i, " x = ", x, " theta_e = ", θe)

        for j = 1:num_states
            y[j]   = solution[i,j];
            ini[j] = solution[i,j];
        end
            ini[num_states+1] = θe;  #Actual θe;

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

        if i==1
            solution[Nspan,9]=solution[Nspan,5]./solution[Nspan,8]
            solution[Nspan,10]=solution[Nspan,3]./solution[Nspan,7]
        end

        #Update vector of states:
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;
        if i < Nspan
            solution[i+1,1:8] = solution[i,1:8] - step*dy;
            solution[i+1,9]=solution[i+1,5]./solution[i+1,7]
            solution[i+1,10]=solution[i+1,3]./solution[i+1,8]
        end
    end
end
