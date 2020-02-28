function marginal_taxes(ctrlvec::Array{Float64},θvec::Array{Float64},solvec::Array{Float64},pa)

    (Nspan,~)=size(solvec)

    for j=1:Nspan
      #Definition of states and controls:
      θ = θvec[j];

      uw    = solvec[j,1];
      μ     = solvec[j,2];
      e     = solvec[j,3];
      ϕ_e   = solvec[j,4];
      y_agg = solvec[j,5];
      λ     = solvec[j,6];
      l_agg = solvec[j,7];
      ω     = solvec[j,8];

      zz    = ctrlvec[j,1];
      nn    = ctrlvec[j,2];
      ll    = ctrlvec[j,3];
      pp    = ctrlvec[j,4];

      #Find the marginal taxes:
      τ_prime[j,1] = pa.β*zz^pa.σ; #τ_c_prime
      τ_prime[j,2] = (λ*pa.α*e*nn^(pa.α-1)-ω)/ω; #τ_n_prime
      τ_prime[j,3] = (θ*ω-λ*pa.χ*ll^(pa.ψ))/(θ*ω); #τ_l_prime

    end
end

function marginal_taxese(ctrlvec::Array{Float64},θvec::Array{Float64},solvec::Array{Float64},pa)

    (Nspan,~)=size(solvec)

    for j=1:Nspan
      #Definition of states and controls:
      θ = θvec[j];

      ue    = solvec[j,1];
      μe    = solvec[j,2];
      ye    = solvec[j,3];
      λe    = solvec[j,4];
      le    = solvec[j,5];
      ωe    = solvec[j,6];

      zz    = ctrlvec[j,1];
      nn    = ctrlvec[j,2];

      #Find the marginal taxes:
      τ_prime_e[j,1] = pa.β*zz^pa.σ; #τ_c_prime
      τ_prime_e[j,2] = (λe*pa.α*θ*nn^(pa.α-1)-ωe)/ωe; #τ_n_prime

    end
end

function taxes_path(ctrlvec::Array{Float64},θvec::Array{Float64},solvec::Array{Float64},pa,θevec::Array{Float64},solvece::Array{Float64},ctrlvece::Array{Float64})

    (Nspan,~)=size(solvec)

    #Defining the taxes in the θ_w_lb: the order: Tc, Tn and Tl (first for θw, then Entrepreneurs θe).
    Taxes    = Array{Float64}(undef,Nspan,5)

    Taxes[1,1] = solvec[1,3]*ctrlvec[1,2]^pa.α-solvec[1,8]/solvec[1,6]*ctrlvec[1,2]-pa.β/(1.0+pa.σ)*ctrlvec[1,1]^(1.0+pa.σ); #T_c
    Taxes[1,2] = 0.0; #T_n
    Taxes[1,3] = θvec[1]*solvec[1,8]/solvec[1,6]*ctrlvec[1,3]-pa.χ/(1.0+pa.ψ)*ctrlvec[1,3]^(1.0+pa.ψ)-solvec[1,1]; #T_l

    #First we solve the taxes for θw:
    for j = 2:Nspan
      #Definition of states and controls:
        #θ:
        θ     = θvec[j];

        uw    = solvec[j,1];
        μ     = solvec[j,2];
        e     = solvec[j,3];
        ϕ_e   = solvec[j,4];
        y_agg = solvec[j,5];
        λ     = solvec[j,6];
        l_agg = solvec[j,7];
        ω     = solvec[j,8];

        zz    = ctrlvec[j,1];
        nn    = ctrlvec[j,2];
        ll    = ctrlvec[j,3];
        pp    = ctrlvec[j,4];

        Taxes[j,1] = pa.β*zz^pa.σ; #Tc
        Taxes[j,2] = (λ*pa.α*e*nn^(pa.α-1.0)-ω)/ω; #Tn
        Taxes[j,3] = (θ*ω-λ*pa.χ*ll^pa.ψ)/(θ*ω); #Tl

    end

    #Now we solve the taxes for θe:
    Taxes[1,4] = Taxes[end,1]; #T_c
    Taxes[1,5] = Taxes[end,2]; #T_n

    for j = 2:Nspan
      #Definition of states and controls:
        #θe:
        θe    = θevec[j];

        ue    = solvece[j,1];
        μe    = solvece[j,2];
        ye    = solvece[j,3];
        λe    = solvece[j,4];
        le    = solvece[j,5];
        ωe    = solvece[j,6];

        zz    = ctrlvece[j,1];
        nn    = ctrlvece[j,2];

        Taxes[j,4] = pa.β*zz^pa.σ; #Tc
        Taxes[j,5] = (λe*pa.α*θe*nn^(pa.α-1.0)-ωe)/ωe; #Tn
    end

    #We have to get the values of the integrals:
    for i=1:5
      to_integrate = Taxes[:,i];
      sol_int = integrals(to_integrate,θspan);
      Taxes[:,i] = sol_int;
    end


    #Final Output:
    Taxes;
end
