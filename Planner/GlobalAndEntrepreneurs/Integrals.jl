function my_integral_lb!(sol_int::Array{Float64,1},intvec::Array{Float64,1},θvec::Array{Float64,1},y0::Float64)

    (Nspan,) = size(θvec);

    #We are going to integrate using the trapezoid rule:
    for i=1:Nspan
        if i == 1
          sol_int[i] = y0;
        else
          sol_int[i] = sol_int[i-1]+(intvec[i]+intvec[i-1])/2.0*(θvec[i]-θvec[i-1]);
        end
    end
    nothing
end

function my_integral_ub!(sol_int::Array{Float64,1},intvec::Array{Float64,1},θvec::Array{Float64,1},yend::Float64)

    (Nspan,) = size(θvec);

    #We are going to integrate using the trapezoid rule:
    for i=Nspan:-1:1
        if i == Nspan
          sol_int[i] = yend;
        else
          sol_int[i] = sol_int[i+1]-(intvec[i]+intvec[i+1])/2.0*(θvec[i]-θvec[i+1]);
        end
    end
    nothing
end
