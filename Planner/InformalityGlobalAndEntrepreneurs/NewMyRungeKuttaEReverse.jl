function my_runge_kuttae_reverse!(solution::Array{Float64,2},y_end::Array{Float64,1},xspan,step::Float64,pa,θw::Float64,controlsRK::Array{Float64,2}; verbose = false)
#The states vector is given in the following order:
# ue, μe, ye, λe, lfe, ωfe, lie, ωie, wie and ϕwe.
# We get the following auxiliary states: le_new, lie_new and ye_new.

    #θw is the upper bound of workers distribution. It´s taken from the global problem.
    (Nspan,columns) = size(solution);
    num_states      = columns - 3;
    solution[end,1:num_states] = y_end;
    solution[end,14] = solution[end,5]./solution[end,11];
    solution[end,15] = solution[end,7]./solution[end,12];
    solution[end,16] = solution[end,3]./solution[end,13];

    #Defining the vectors for Ruge-Kutta:
    z1  = Array{Float64,1}(undef,num_states);
    z2  = Array{Float64,1}(undef,num_states);
    z3  = Array{Float64,1}(undef,num_states);
    z4  = Array{Float64,1}(undef,num_states);

    for i ∈ [z1,z2,z3,z4]
        fill!(i,NaN);
    end

    #Loop over values of θ_e:
    for i=Nspan:-1:1

        println("i = ", i)
        #Current value for θ_e
        θe = xspan[i];
        println(θe)

        #We save the approximation of the derivative of the states, according to the order 4 Runge-Kutta method:
        #First we need to calculate the optimal controls we are using to solve the states:
        sse = StateE(solution[i,1], solution[i,2], solution[i,3], solution[i,4], solution[i,5], solution[i,6], solution[i,7], solution[i,8], solution[i,9], solution[i,10]); #State object.

        (ze, ne, nie) = new_find_controlse(θe, sse, pa);
        println("ze = ", ze, "ne = ", ne, "nie = ", nie)
        any(isnan,(ze, ne, nie)) && error("Function find_statese gets NaN controls.")
        controls = [ze, ne, nie];

        #println("z1")
        find_statese!(z1, solution[i,1:num_states], pa, θe, controls);
        any(isnan,z1) && error("z1 is NaN ")
        #println("z2")
        find_statese!(z2, solution[i,1:num_states].+0.5*step*z1, pa, θe+0.5*step, controls);
        any(isnan,z2) && error("z2 is NaN ")
        #println("z3")
        find_statese!(z3, solution[i,1:num_states].+0.5*step*z2, pa, θe+0.5*step, controls);
        any(isnan,z3) && error("z3 is NaN ")
        #println("z4")
        find_statese!(z4, solution[i,1:num_states].+step*z3, pa, θe+step, controls);
        any(isnan,z4) && error("z4 is NaN ")

        #Update vector of states
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;
        if i > 1
            solution[i-1,1:num_states] = solution[i,1:num_states] - step*dy;

            solution[i-1,14] = solution[i-1,5]./solution[i-1,11];
            solution[i-1,15] = solution[i-1,7]./solution[i-1,12];
            solution[i-1,16] = solution[i-1,3]./solution[i-1,13];
        end
    end
end
