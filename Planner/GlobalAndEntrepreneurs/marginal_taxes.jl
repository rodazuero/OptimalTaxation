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

function taxes_path(ctrlvec::Array{Float64},θvec::Array{Float64},solvec::Array{Float64},pa)

    (Nspan,~)=size(solvec)

    #Defining the taxes in the θ_w_lb:
    TaxesTc    = Array{Float64}(undef,Nspan,1)
    TaxesTn    = Array{Float64}(undef,Nspan,1)
    TaxesTl    = Array{Float64}(undef,Nspan,1)

    TaxesTc[1] = solvec[1,3]*ctrlvec[1,2]^pa.α-solvec[1,8]/solvec[1,6]*ctrlvec[1,2]-pa.β/(1.0+pa.σ)*ctrlvec[1,1]^(1.0+pa.σ); #T_c
    TaxesTn[1] = 0.0; #T_n
    TaxesTl[1] = θvec[1]*solvec[1,8]/solvec[1,6]*ctrlvec[1,3]-pa.χ/(1.0+pa.ψ)*ctrlvec[1,3]^(1.0+pa.ψ)-solvec[1,1]; #T_l

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

        TaxesTc[j] = pa.β*zz^pa.σ;
        TaxesTn[j] = (λ*pa.α*e*nn^(pa.α-1.0)-ω)/ω;
        TaxesTl[j] = (θ*ω-λ*pa.χ*ll^pa.ψ)/(θ*ω);

    end

    Int_tax      = Array{Float64}(undef,3,1);
    Int_tax[1,1] = (sum(TaxesTc)-0.5*TaxesTc[1]-0.5*TaxesTc[end])*(θvec[end]-θvec[1])/(Nspan-1.0)*(1.0-pa.constant_w_lw*pa.constant_e_lw);
    Int_tax[2,1] = (sum(TaxesTn)-0.5*TaxesTn[1]-0.5*TaxesTn[end])*(θvec[end]-θvec[1])/(Nspan-1.0)*(1.0-pa.constant_w_lw*pa.constant_e_lw);
    Int_tax[3,1] = (sum(TaxesTl)-0.5*TaxesTl[1]-0.5*TaxesTl[end])*(θvec[end]-θvec[1])/(Nspan-1.0)*(1.0-pa.constant_w_lw*pa.constant_e_lw);

    #Final Output:
    Int_tax, TaxesTc, TaxesTn, TaxesTl;
end

function taxes_pathe(ctrlvec::Array{Float64},θvec::Array{Float64},solvec::Array{Float64},pa)

    (Nspan,~)=size(solvec)

    #Defining the elasticities and the taxes in the θ_w_lb:
    h = θvec[2]-θvec[1]; #Distance between each element.

    TaxesE[1,1]  = solvec[1,3]*ctrlvec[1,2]^pa.α-solvec[1,8]/solvec[1,6]*ctrlvec[1,2]-pa.β/(1.0+pa.σ)*ctrlvec[1,1]^(1.0+pa.σ); #T_c
    TaxesE[1,2]  = 0.0; #T_n

    for j=1:Nspan-1
      #Definition of states and controls:
        #θ:
        θ     = θvec[j];

        ue    = solvec[j,1];
        μe    = solvec[j,2];
        ye    = solvec[j,3];
        λe    = solvec[j,4];
        le    = solvec[j,5];
        ωe    = solvec[j,6];

        zz    = ctrlvec[j,1];
        nn    = ctrlvec[j,2];

        #θ+1:
        θ1 = θvec[j+1];

        ue1    = solvec[j+1,1];
        μe1    = solvec[j+1,2];
        ye1    = solvec[j+1,3];
        λe1    = solvec[j+1,4];
        le1    = solvec[j+1,5];
        ωe1    = solvec[j+1,6];

        zz1    = ctrlvec[j+1,1];
        nn1    = ctrlvec[j+1,2];

        TaxesE[j+1,1]  = TaxesE[j,1]+h/2.0*(pa.β*zz^pa.σ + pa.β*zz1^pa.σ);
        TaxesE[j+1,2]  = TaxesE[j,2]+h/2.0*((λe*pa.α*θ*nn^(pa.α-1.0)-ωe)/ωe + (λe1*pa.α*θ1*nn1^(pa.α-1.0)-ωe1)/ωe1);

    end

end
