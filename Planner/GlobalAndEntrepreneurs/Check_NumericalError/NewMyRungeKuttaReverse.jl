
function my_runge_kutta_reverse!(solution::Array{Float64},y_end,xspan,step, pa, par, cons; verbose = false)
    #
    solution[end,:] = y_end

    z1 = Array{Float64,1}(undef,2)
    z2 = Array{Float64,1}(undef,2)
    z3 = Array{Float64,1}(undef,2)
    z4 = Array{Float64,1}(undef,2)

    y = Array{Float64,1}(undef,2)
    ini = Array{Float64,1}(undef,3)

    (Nspan,) =  size(xspan)

    #Loop over values of θw
    for i = 1:Nspan

        println("i = ", i)
        #Current value for theta
        x = xspan[i];
        θ = xspan[i];

        #Convert state object into a vector
        y[1] = solution[Nspan+1-i,1]  # μ
        y[2] = solution[Nspan+1-i,2]  # e

        ini[1] = solution[Nspan+1-i,1]  # μ
        ini[2] = solution[Nspan+1-i,2]  # e
        ini[3] = θ  #Actual θ_w;

        find_states!(z1, y, pa, θ, ini, par, cons);
        any(isnan,z1) && error("z1 is NaN ")

        #println("z2")
        find_states!(z2, y+0.5*step*z1, pa, x+0.5*step, ini, par, cons);
        any(isnan,z2) && error("z2 is NaN ")

        #println("z3")
        find_states!(z3, y+0.5*step*z2, pa, x+0.5*step, ini, par, cons);
        any(isnan,z3) && error("z3 is NaN ")

        #println("z4")
        find_states!(z4, y+step*z3, pa, x+step, ini, par, cons);
        any(isnan,z4) && error("z4 is NaN ")

        #Update vector of states
        dy = (z1+2.0*z2+2.0*z3+z4)/6.0;
        if i < Nspan
            solution[Nspan-i,1:2] = solution[Nspan+1-i,1:2] - step*dy;
        end

    end

end
