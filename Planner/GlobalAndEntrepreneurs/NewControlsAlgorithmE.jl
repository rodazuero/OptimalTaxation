function new_find_controlse!( controlse::Array{Float64,1}, θe::Float64, sse::StateE, pa, verbose )
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    # 0.0 Preallocating
        current_z::Float64 = 0.0;
        current_n::Float64 = 0.0;
        current_ham::Float64 = 0.0;
        candidate_z::Float64 = NaN;
        candidate_n::Float64 = NaN;
        candidate_ham::Float64 = NaN;
        previousCPOpositive::Bool = false; #Indicates change from positive to negative in the FOC.
        candidate_cpo::Float64 = NaN;

    # 0.1 Recover densities:
        h_e= pa.he(pa.θ_w_ub, θe);

    # 0.2 Initialize constants
        n_full_info::Float64 = ((sse.λe*pa.α*θe)/sse.ωe)^(1.0/(1.0-pa.α));
        z_max::Float64 = (1.0/pa.β)^(1.0/pa.σ);
        n_lwbar::Float64 = pa.ς; #The smallest possible n.
        n_upbar::Float64 = n_full_info;

    # 0.3 Define auxiliar functions:
        z_opt(nvar)     = - pa.σ*sse.μe*nvar^pa.α/(sse.λe*h_e);
        nfoc(nvar,zvar) = sse.λe*pa.α*θe*nvar^(pa.α-1.0)*h_e - sse.ωe*h_e + pa.α*sse.μe*nvar^(pa.α-1.0)*(1.0-pa.β*zvar^pa.σ);
        nfoc_at_zopt(nvar) = nfoc(nvar,z_opt(nvar))
        objective(nvar,zvar) = ( pa.indicator*sse.ue^pa.ϕ*h_e +
                                 sse.λe*( θe*nvar^pa.α - pa.β*zvar^(1.0+pa.σ)/(1.0+pa.σ) - sse.ue)*h_e -
                                 sse.ωe*nvar*h_e + sse.μe*nvar^pa.α*(1.0-pa.β*zvar^pa.σ) );

     # 0.4 For Grid
         Nz::Integer = 100; #Number of Ns
         nstep::Float64 = (n_upbar - n_lwbar)/(Nz-1);
         verbose && println("nstep = ", nstep, ", Nz = ", Nz)
         ngrid = range(n_lwbar, stop = n_upbar, length = Nz);

    # 1. Define the grid for the algorithm, we don't have a monotone function, so we cannot run the bisection directly:
    for j = 1:Nz

        candidate_cpo = nfoc_at_zopt(ngrid[j]);

        if previousCPOpositive == 1
            if candidate_cpo < 0
                candidate_n   = find_zero(nfoc_at_zopt, (ngrid[j-1],ngrid[j]), Bisection());
                candidate_z   = z_opt(candidate_n);
                candidate_ham = objective(candidate_n,candidate_z);
                    if candidate_ham > current_ham
                        current_n   = candidate_n;
                        current_z   = candidate_z;
                        current_ham = candidate_ham;
                    end #if
                previousCPOpositive = false;
            end #if
        else
                candidate_cpo > 0 ? previousCPOpositive = true : previousCPOpositive = false;
        end #if
     end #for

     if candidate_cpo>0 #here the value corresponds to the last value in the grid (n_upbound)
         candidate_n   = n_upbar;
         candidate_z   = z_opt(candidate_n);
         candidate_ham = objective(candidate_n,candidate_z);
         if candidate_ham > current_ham
             current_n  = candidate_n;
             current_z  = candidate_z;
             current_ham = candidate_ham;
         end #if
     end #if

    #Final Output:
    controlse[1] = current_n;
    controlse[2] = current_z;
    #println("zze = ", zze, "nne = ", nne)

end

function new_find_controlse( θe::Float64, sse::StateE, pa, verbosebool::Bool = false)
    controlse = Array{Float64}(undef,2);
    new_find_controlse!(controlse, θe::Float64, sse::StateE, pa, verbosebool);
    return Tuple(controlse)
end

function recover_controlse!(ctrlvec::Array{Float64}, θvec::Array{Float64}, solvec::Array{Float64})
    (Nspan,~)=size(solvec)

    for j = Nspan:-1:1
      θe  = θvec[j];
      ue  = solvec[j,1];
      μe  = solvec[j,2];
      λe  = solvec[j,4];
      ωe  = solvec[j,6];
      sse = StateE(ue, μe, λe, ωe);

      (nne, zze)   = new_find_controlse( θe, sse, pa)
      ctrlvec[j,1] = zze;
      ctrlvec[j,2] = nne;
    end
    nothing
end
