function integrals(intvec::Array{Float64},θvec::Array{Float64},method)
#Method indicates what limit in the integral we have fixed:
    #1: we have the upper limit fixed,
    #0: we have the lower limit fixed.

    (Nspan,) = size(θvec);
    sol_int  = fill(NaN,Nspan,1);
    factor   = (θvec[end]-θvec[1])/(Nspan-1.0); #Factor to multiply the integral.

    #We are going to integrate using the trapezoid rule:
    if method == 1
        #When we have the upper bound fixed in the integral:
        for i=1:Nspan
            sol_int[i,1] = (sum(intvec[i:end])-0.5*intvec[i]-0.5*intvec[end])*factor;
        end
    else
        #When we have the lowed bound fixed in the integral:
        for i=1:Nspan
            sol_int[i,1] = (sum(intvec[1:i])-0.5*intvec[i]-0.5*intvec[1])*factor;
        end
    end

    #Returns of the function --> Vector containing the integral for each productivity:
    sol_int;
end
