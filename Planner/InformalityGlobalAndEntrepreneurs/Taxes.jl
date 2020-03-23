function marginal_taxes!(τ_prime::Array{Float64,2},solvec::Array{Float64,2},ctrlvec::Array{Float64,2},θvec::Array{Float64,1},pa)

    (Nspan,~)=size(solvec)

    for j=Nspan:-1:1
      #Definition of states and controls:
      θ = θvec[j];

      uw     = solvec[1];
      μ      = solvec[2];
      e      = solvec[3];
      ϕe     = solvec[4];
      y_agg  = solvec[5];
      λ      = solvec[6];
      lf_agg = solvec[7];
      ωf     = solvec[8];
      li_agg = solvec[9];
      ωi     = solvec[10];
      wi     = solvec[11];
      ϕw     = solvec[12];

      zz    = ctrlvec[j,1];
      nn    = ctrlvec[j,2];
      nni   = ctrlvec[j,3];
      ll    = ctrlvec[j,4];
      lli   = ctrlvec[j,5];
      pp    = ctrlvec[j,6];

      #Find the marginal taxes:
      τ_prime[j,1] = pa.β*zz^pa.σ; #τc prime
      τ_prime[j,2] = (λ*pa.δ*nni^pa.γ+λ*wi-ωf)/ωf; #τn prime
      τ_prime[j,3] = (ωf - θ*wi - pa.κ*θ^(1.0+pa.ρ)*li^pa.ρ)/ωf; #τl prime
    end
end

function marginal_taxese!(τ_prime_e::Array{Float64,2},solvec::Array{Float64,2},ctrlvec::Array{Float64,2},θvec::Array{Float64,1},pa)

    (Nspan,~)=size(solvec)

    for j=Nspan:-1:1
      #Definition of states and controls:
      θe = θvec[j];

      ue  = solvec[1];
      μe  = solvec[2];
      ye  = solvec[3];
      λe  = solvec[4];
      lfe = solvec[5];
      ωfe = solvec[6];
      lie = solvec[7];
      ωie = solvec[8];
      wie = solvec[9];
      ϕwe = solvec[10];

      zze = ctrlvec[j,1];
      nne = ctrlvec[j,2];
      nni = ctrlvec[j,2];

      #Find the marginal taxes:
      τ_prime_e[j,1] = pa.β*zze^pa.σ; #τc prime
      τ_prime_e[j,2] = (λe*θe*nne^(pa.α-1.0) - ωfe)/ωfe; #τn prime
      τ_prime_e[j,3] = (λe*pa.δ*nni^pa.γ+λe*wie-ωfe)/ωfe; #τn prime (formula from ni)

    end
end








function taxes_path(ctrlvec::Array{Float64},θvec::Array{Float64},solvec::Array{Float64},pa,θevec::Array{Float64},solvece::Array{Float64},ctrlvece::Array{Float64},margtax::Array{Float64},margtaxe::Array{Float64})

    (Nspan,~)=size(solvec)

    #Defining the taxes in the θ_w_lb: the order: Tc, Tn and Tl (first for θw, then Entrepreneurs θe).
    #Taxes    = Array{Float64}(undef,Nspan,5)
    Taxes_gb      = fill(NaN,Nspan,3)
    Taxes_gb[1,1] = solvec[1,3]*ctrlvec[1,2]^pa.α-solvec[1,8]/solvec[1,6]*ctrlvec[1,2]-pa.β/(1.0+pa.σ)*ctrlvec[1,1]^(1.0+pa.σ)-solvec[1,1]; #T_c
    Taxes_gb[1,2] = 0.0; #T_n
    Taxes_gb[1,3] = θvec[1]*solvec[1,8]/solvec[1,6]*ctrlvec[1,3]-pa.χ/(1.0+pa.ψ)*ctrlvec[1,3]^(1.0+pa.ψ)-solvec[1,1]; #T_l

    #We have to get the values of the integrals:
    #Define the values of the limits of the integral (tax base):
    tax_base = fill(NaN,Nspan,3)
    for j = 1:Nspan
      θ  = θspan[j];
      e  = solvec[j,3];
      ω  = solvec[j,8];
      nn = ctrlvec[j,2];
      ll = ctrlvec[j,3];
      tax_base[j,2] = ω*nn; #Tn
      tax_base[j,3] = θ*ll*ω; #Tl
    end

    for i=2:3
      sol_int      = fill(NaN,Nspan)
      base         = tax_base[:,i]; #Tax base for the integral factor
      to_integrate = margtax[:,i];
      my_integral_lb!(sol_int,to_integrate,base,Taxes_gb[1,i]);
      Taxes_gb[:,i] = sol_int[:];
    end

    #Calculating for Tc:
    for j = 1:Nspan
      θ  = θspan[j];
      e  = solvec[j,3];
      ω  = solvec[j,8];
      zz = ctrlvec[j,1];
      nn = ctrlvec[j,2];
      ll = ctrlvec[j,3];
      Tn = Taxes_gb[j,2];
      tax_base[j,1] = e*nn^pa.α-ω*nn-Tn-zz; #Tc
    end

    i=1
      sol_int      = fill(NaN,Nspan)
      base         = tax_base[:,i]; #Tax base for the integral factor
      to_integrate = margtax[:,i];
      my_integral_lb!(sol_int,to_integrate,base,Taxes_gb[1,i]);
      Taxes_gb[:,i] = sol_int[:];

    #Now we solve the taxes for θe:
    Taxes_ent      = fill(NaN,Nspan,2)
    Taxes_ent[1,1] = Taxes_gb[end,1]; #T_c
    Taxes_ent[1,2] = Taxes_gb[end,2]; #T_n

    tax_base_ent = fill(NaN,Nspan,2)
    for j = 1:Nspan
      ωe = solvece[j,6];
      nn = ctrlvece[j,2];
      tax_base_ent[j,2] = ωe*nn; #Tn
    end

    i=2
      sol_int      = fill(NaN,Nspan)
      base         = tax_base_ent[:,i]; #Tax base for the integral factor
      to_integrate = margtaxe[:,i];
      my_integral_lb!(sol_int,to_integrate,base,Taxes_ent[1,i]);
      Taxes_ent[:,i] = sol_int[:];

    #Calculating for Tc:
    for j = 1:Nspan
      θe = θespan[j];
      ωe = solvece[j,6];
      zz = ctrlvece[j,1];
      nn = ctrlvece[j,2];
      Tn = Taxes_ent[j,2]
      tax_base_ent[j,1] = θe*nn^pa.α-ωe*nn-Tn-zz; #Tc
    end

    i=1
      sol_int      = fill(NaN,Nspan)
      base         = tax_base_ent[:,i]; #Tax base for the integral factor
      to_integrate = margtaxe[:,i];
      my_integral_lb!(sol_int,to_integrate,base,Taxes_ent[1,i]);
      Taxes_ent[:,i] = sol_int[:];

    #Final Output:
    Taxes_gb, Taxes_ent, tax_base, tax_base_ent;
end
