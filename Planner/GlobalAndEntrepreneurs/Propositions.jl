function propositions(ctrlvec::Array{Float64},θvec::Array{Float64},solvec::Array{Float64},tax_prime::Array{Float64},pa,θevec::Array{Float64},solvece::Array{Float64})

    (Nspan,~) = size(solvec)
    (Nespan,) = size(θevec)

    #Defining the propostions matrix:
    proposition1 = Array{Float64}(undef,Nspan,3)
    proposition2 = Array{Float64}(undef,Nspan,2)
    proposition3 = Array{Float64}(undef,Nspan,4)

    #Defining the elasticities and the taxes in the θ_w_lb:
    ε_l_1Tlpr = 1.0/pa.ψ;
    ε_l_θw    = 1.0/pa.ψ;
    ε_z_Tcpr  = 1.0/pa.σ;

    T_c = tax_prime[1]
    T_n = tax_prime[2]
    T_l = tax_prime[3]

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
      h_e     = pa.he(θ, e);
      h_e_fix = pa.he(pa.θ_w_ub, e);
      h_w     = pa.hw(θ, e);
      hh      = pa.hh(θ, e, pp);
      gg      = pa.gg(θ, e);

      T_cpr = tax_prime[j,1];
      T_npr = tax_prime[j,2];
      T_lpr = tax_prime[j,3];

      Vw = (pa.indicator*uw - λ*uw) + ω*ll*θ - λ*pa.χ/(1.0+pa.ψ)*ll^(1.0+pa.ψ);
      Ve = (pa.indicator*uw - λ*uw) + λ*e*nn^pa.α - λ*pa.β/(1.0+pa.σ)*zz^(1.0+pa.σ) - ω*nn;

      #Solving the elements of the propositions:
      #1:
        proposition1[j,1] = (1.0-pa.indicator*pa.ϕ*uw^(pa.ϕ-1.0)/λ)*hh;
        proposition1[j,2] = T_lpr/(1.0-T_lpr)*ε_l_1Tlpr/(1.0+ε_l_1Tlpr)*θ*h_w;
        proposition1[j,3] = ε_z_Tcpr*zz/nn^pa.α*h_e;

        #2:
        proposition2[j,1] = ε_z_Tcpr*zz/nn^pa.α;
        proposition2[j,2] = e/(1.0-T_cpr)*T_npr/(1.0+T_npr);

        #3:
        proposition3[j,1] = ε_z_Tcpr*zz/nn^pa.α*h_e;
        proposition3[j,2] = 1.0/λ*(Vw-Ve)*gg*1.0/(nn^pa.α*(1.0-pa.β*zz^pa.σ));
        proposition3[j,3] = (1.0-pa.indicator*pa.ϕ*uw^(pa.ϕ-1.0)/λ)*h_e_fix;

    end

    for j=1:Nespan
      #Definition of states and controls:
      θe = θevec[j];

      ue    = solvece[j,1];
      μe    = solvece[j,2];
      ye    = solvece[j,3];
      λe    = solvece[j,4];
      le    = solvece[j,5];
      ωe    = solvece[j,6];

      #Recover ditributions and elasticities:
      h_e_fix = pa.he(pa.θ_w_ub, θe);

        #3:
        proposition3[j,4] = (1.0-pa.indicator*pa.ϕ*ue^(pa.ϕ-1.0)/λe)*h_e_fix;

    end

    #Returns of the function:
    proposition1, proposition2, proposition3;

end
