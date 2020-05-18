function new_find_controlse!( controlse::Array{Float64,1}, θe::Float64, sse::StateE, pa, verbose)
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    # 0.0 Preallocating
        corner_n::Float64 = NaN
    	corner_z::Float64 = NaN
    	corner_hamiltonian::Float64 = NaN
        current_hamiltonian::Float64 = NaN

    # 0.1 Recover densities:
        h_e= pa.he(pa.θ_w_ub,θe);

    # 0.2 Initialize constants
        n_full_info::Float64 = ((sse.λe*pa.α*θe)/sse.ωe)^(1.0/(1.0-pa.α));
        z_max::Float64 = (1.0/pa.β)^(1.0/pa.σ);
        n_lwbar::Float64 = pa.ς; #The smallest possible n.
        n_upbar::Float64 = ( - sse.λe*z_max*h_e / (sse.μe*pa.σ))^(1.0/pa.α);
        n_upbar = min(n_full_info,n_upbar);

    # 0.3 Define auxiliar functions:
        z_opt(nvar)     = - pa.σ*sse.μe*nvar^pa.α/(sse.λe*h_e);
        nfoc(nvar,zvar) = sse.λe*pa.α*θe*nvar^(pa.α-1.0)*h_e - sse.ωe*h_e + pa.α*sse.μe*nvar^(pa.α-1.0)*(1.0-pa.β*zvar^pa.σ);
        nfoc_at_zopt(nvar) = nfoc(nvar,z_opt(nvar))
        objective(nvar,zvar) = ( pa.indicator*sse.ue^pa.ϕ*h_e +
                                 sse.λe*( θe*nvar^pa.α - pa.β*zvar^(1.0+pa.σ)/(1.0+pa.σ) - sse.ue)*h_e -
                                 sse.ωe*nvar*h_e + sse.μe*nvar^pa.α*(1.0-pa.β*zvar^pa.σ) );

    # 1. Define the algorithm:
    if sse.μe < 1e-15
        #In θe = θ_e_ub μe = 0.0, then we can solve analytically the problem:
        controlse[2] = 0.0;
        controlse[1] = n_full_info;
    else
        #For other values we can solve using bisection:
        if nfoc_at_zopt(n_lwbar)*nfoc_at_zopt(n_upbar)<0
            controlse[1] = find_zero(nfoc_at_zopt, (n_lwbar,n_upbar), Bisection());
            controlse[2] = z_opt(controlse[1]);
            current_hamiltonian = objective(controlse[1],controlse[2]);
        end #if

        #Now we get the corner solution (if any): First when n = ς, and then when n = n_full_info
        corner_n = n_lwbar;
        corner_z = z_opt(corner_n);
        corner_hamiltonian = objective(corner_n,corner_z);
        if corner_hamiltonian > current_hamiltonian
            controlse[1]   = corner_n;
            controlse[2]   = corner_z;
            current_hamiltonian = corner_hamiltonian;
        end #if

        corner_n = n_upbar;
        corner_z = z_opt(corner_n);
        corner_ham = objective(corner_n,corner_z);
        if corner_hamiltonian > current_hamiltonian
            controlse[1]   = corner_n;
            controlse[2]   = corner_z;
            current_hamiltonian = corner_hamiltonian;
        end #if
    end #if

    #Final Output:
    verbose && println("zze = ", controlse[2], "nne = ", controlse[1])
end

function new_find_controlse( θe::Float64, sse::StateE, pa, verbosebool::Bool = false)
    controlse = Array{Float64}(undef,2)
    new_find_controlse!(controlse, θe::Float64, sse::StateE, pa, verbosebool)
    return Tuple(controlse)
end

function recover_controlse!(ctrlvec::Array{Float64},θvec::Array{Float64}, solvec::Array{Float64})
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
