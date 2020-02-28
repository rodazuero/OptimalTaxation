function integrals(intvec::Array{Float64},θvec::Array{Float64})

    (Nspan,) = size(θvec);
    sol_int = Array{Float64}(undef,Nspan,1);
    
    #We are going to integrate using the trapezoid rule:
    for i=1:Nspan
        sol_int[i,1] = (sum(intvec[i:Nspan])-0.5*intvec[i]-0.5*intvec[end])*(θvec[end]-θvec[1])/(Nspan-1.0);
    end

    #Returns of the function --> Vector containing the integral for each productivity:
    sol_int;
end
