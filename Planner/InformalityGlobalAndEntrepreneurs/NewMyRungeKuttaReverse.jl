function my_runge_kuttae_reverse!(solution::Array{Float64,2},y_end::Array{Float64,1},xspan,step::Float64,pa,θw::Float64,controlsRK::Array{Float64,2}; verbose = false)
#The states vector is given in the following order:
# uw, μ, e, ϕe, y, λ, lf, ωf, li, ωi, wi and ϕw.
# We get the following auxiliary states: lf_new, li_new and y_new.

    #θw is the upper bound of workers distribution. It´s taken from the global problem.
    (Nspan,columns) = size(solution);
    num_states      = columns - 3;
    solution[end,1:num_states] = y_end;
    solution[end,16] = solution[end,5]./solution[end,11]; #Lf/Lf_accumualted
    solution[end,17] = solution[end,7]./solution[end,12]; #Li/Li_accumualted
    solution[end,18] = solution[end,3]./solution[end,13]; #Y/Y_accumualted

    #Defining the vectors for Ruge-Kutta:
    z1  = Array{Float64,1}(undef,num_states);
    z2  = Array{Float64,1}(undef,num_states);
    z3  = Array{Float64,1}(undef,num_states);
    z4  = Array{Float64,1}(undef,num_states);

    for i ∈ [z1,z2,z3,z4]
        fill!(i,NaN);
    end

    #Loop over values of θw:
    for i=Nspan:-1:1

        println("i = ", i)
        #Current value for θw
        θe = xspan[i];
        println(θw)

        #We save the approximation of the derivative of the states, according to the order 4 Runge-Kutta method:
        #First we need to calculate the optimal controls we are using to solve the states:
        sse = State(solution[i,1], solution[i,2], solution[i,3], solution[i,4], solution[i,5], solution[i,6], solution[i,7], solution[i,8], solution[i,9], solution[i,10], solution[i,11], solution[i,12]); #State object.

        (z, n, ni, l, li, p) = new_find_controls(θ, sse, pa);
        println("z = ", z, "n = ", n, "ni = ", ni, "l = ", l, "li = ", li)
        any(isnan,(z, n, ni, l, li)) && error("Function find_statese gets NaN controls.")
        controls = [z, n, ni, l, li, p];

        #println("z1")
        find_states!(z1, solution[i,1:num_states], pa, θe, controls);
        any(isnan,z1) && error("z1 is NaN ")
        #println("z2")
        find_states!(z2, solution[i,1:num_states].+0.5*step*z1, pa, θe+0.5*step, controls);
        any(isnan,z2) && error("z2 is NaN ")
        #println("z3")
        find_states!(z3, solution[i,1:num_states].+0.5*step*z2, pa, θe+0.5*step, controls);
        any(isnan,z3) && error("z3 is NaN ")
        #println("z4")
        find_states!(z4, solution[i,1:num_states].+step*z3, pa, θe+step, controls);
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
