
function my_runge_kutta_reverse!(solution::Array{Float64},y_end,xspan,step,pa; verbose = false)

    # θw is the upper bound of workers distribution. It´s taken from the global problem.
    solution[end,1:10] = y_end[1:10]
    solution[end,11]  = solution[end,5]/solution[end,8];
    solution[end,12] = solution[end,3]/solution[end,7];

    #Allocate memory. I follow the notation in Juddm page 345.
    z1 = Array{Float64,1}(undef,10);
    z2 = Array{Float64,1}(undef,10);
    z3 = Array{Float64,1}(undef,10);
    z4 = Array{Float64,1}(undef,10);

    y   = Array{Float64,1}(undef,10);
    ini = Array{Float64,1}(undef,11);

    (Nspan,) =  size(xspan);

    #Loop over values of θ_w:
    for i=Nspan:-1:1

        println("i = ", i)
        #Current value for θ_w
        x = xspan[i];
        θ = xspan[i]
        #println("it = ", i, " x = ", x, " theta_e = ", θe)

        #Convert state object into a vector
        y[:]   = solution[i,1:10]
        ini[1:10] = solution[i,1:10]

        ini[end] = θ  #Actual θw;
        println("θwRK = ", ini[end])

        #We save the approximation of the derivative of the states, according to the order 4 Runge-Kutta method
        #println("z1")
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

        #Update vector of states
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;
        if i > 1
            solution[i-1,1:10] = solution[i,1:10] - step*dy;
            solution[i-1,11]  = solution[i-1,7]./solution[i-1,9]
            solution[i-1,12] = solution[i-1,5]./solution[i-1,10]
        end
    end
end
