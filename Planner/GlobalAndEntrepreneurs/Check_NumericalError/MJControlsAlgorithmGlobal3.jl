function new_find_controls( θ, ss, pa, par, cons)
      
      pp = (1.0+par)*cons*(θ-pa.θ_w_lb)^par

      pp;
end

function recover_controls!(ctrlvec::Array{Float64}, θvec::Array{Float64}, solvec::Array{Float64}, par, cons)
    (Nspan,~)=size(solvec)

    for j=1:Nspan
      θ  = θvec[j];
      μ  = solvec[j,1];
      e  = solvec[j,2];
      #ss = State(μ, e);
      ss = State(e, 0.0, 0.0, μ, 0.0, 0.0)

      #(zz, nn, ll, pp) = new_find_controls( θ, ss, pa)
      (pp) = new_find_controls( θ, ss, pa, par, cons)
      ctrlvec[j,1] = pp;
    end
    nothing
end
