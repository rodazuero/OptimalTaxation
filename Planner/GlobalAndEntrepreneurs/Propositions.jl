function propositions(ctrlvec::Array{Float64},θvec::Array{Float64},solvec::Array{Float64},τ_prime::Array{Float64},pa)

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

      #Recover ditributions and elasticities:
      h_e = pa.he(θ, e);
      h_w = pa.hw(θ, e);

      T_lpr = τ_prime[i,3]
      el_z_tc = 1.0/pa.σ;
      el_l_tl = 1.0/pa.ψ;
      el_l_θw = 1.0/pa.ψ;

      #Solving the elements of the propositions:
      #1:
      T_lpr/(1.0-T_lpr)*el_l_tl/(1.0+el_l_θw)*θ*h_w
      #2:
      e/(1.0-T_cpr)*T_npr/(1.0-T_npr)

      #3:
      el_z_tc*zz/nn^pa.α*h_e

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
