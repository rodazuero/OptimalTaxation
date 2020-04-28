function new_find_controls!(controls::Array{Float64,1}, θ::Float64, ss::State, pa, verbose=false)
	
    #INPUT: states and parameters
    #OUTPUT: optimal controls
# 0.0 Preallocating
    corner_nn::Float64=NaN
	corner_zz::Float64=NaN
	corner_ll::Float64=NaN
	corner_pp::Float64=NaN
	corner_hamiltonian::Float64=NaN
	max_hamiltonian::Float64=-Inf
    verbose && println("θ: ", θ,", e: ",ss.e )
# 0.1 Recover densities:
    h_e::Float64 = pa.he(θ, ss.e);
    h_w::Float64 = pa.hw(θ, ss.e)

# 0.2 Initialize constants
    n_full_info::Float64 = ((ss.λ*pa.α*ss.e)/ss.ω)^(1.0/(1.0-pa.α));
    z_lwbar::Float64 	= 0.0
    z_upbar::Float64 	= (1.0/pa.β)^(1.0/pa.σ)*(1 - 1e-12) #Max possible evasion.
    A_cons::Float64 	= ss.ω*pa.ς+ pa.indicator*ss.uw^pa.ϕ-ss.λ*ss.uw+ss.ϕ_e/h_e;
    nn_z0::Float64 		= -A_cons*1.0/ss.ω*pa.α/(1.0-pa.α);
	pre_tax_profits_at_nmin::Float64 = ss.λ*ss.e*pa.ς^pa.α # - ss.ω*pa.ς

# 0.3 Define aux functions:
    n_opt(zvar)			= pa.α/ss.ω/(1.0-pa.α)*(ss.λ*pa.β/(1.0+pa.σ)*zvar^(1.0+pa.σ)-A_cons)
    zfoc(nvar, zvar)	= ss.λ*ss.e*nvar^pa.α - ss.ω*nvar - ss.λ/pa.σ*zvar*(1.0-pa.β/(1.0+pa.σ)*zvar^pa.σ) + A_cons
    zfoc_at_nopt(zvar)	= zfoc(n_opt(zvar), zvar)
	zfoc_at_nmin(zvar) 	= zfoc(pa.ς, zvar)
    den_l(nvar,zvar)	= ( ss.λ*pa.χ*h_w - pa.χ/θ*(1.0+pa.ψ)*( ss.μ + h_e/( nvar^pa.α*(1.0-pa.β*zvar^pa.σ) )*
    						( ss.λ*ss.e*nvar^pa.α - ss.λ*pa.β/(1.0+pa.σ)*zvar^(1.0+pa.σ) -ss.ω*nvar + A_cons) ) )
    ll_opt(nvar,zvar) 	= (ss.ω*θ*h_w/den_l(nvar,zvar))^(1.0/pa.ψ)
    den_p(nvar,zvar)  	= θ*nvar^pa.α*(1.0-pa.β*zvar^pa.σ)

    objective(lvar, nvar, pvar, zvar) = ( pa.indicator*ss.uw^pa.ϕ*(h_w+pvar*h_e) + ss.μ*pa.χ/θ*lvar^(1.0+pa.ψ)
						    	+ ss.λ*pvar*h_e*( ss.e*nvar^pa.α - pa.β/(1.0+pa.σ)*zvar^(1.0+pa.σ) - ss.uw )
    							- ss.λ*h_w*( ss.uw + pa.χ/(1.0+pa.ψ)*lvar^(1.0+pa.ψ) )
    							+ ss.ω*(θ*lvar*h_w-(nvar-pa.ς)*pvar*h_e) + ss.ϕ_e*pvar ) # Hamiltonian. big parenthesis needed for multline

# 1. Define cases
    # verbose && println("Inside the loop")
    if nn_z0 >= n_full_info
    # Case 0: n_opt(z=0)>n_full_info.
        verbose && println("Case 0: n_opt(z=0)>n_full_info.")
        controls[4] = 0.0 # zz
        controls[2] = nn_z0 # nn
        den_l(controls[2],controls[4]) <= 0.0 && error("Denominator of L is negative or zero.") #Stop when denominator is negative.
        controls[1] = ll_opt(controls[2],controls[4]) 
        controls[3] = (pa.χ*controls[1]^(1.0+pa.ψ))/den_p(controls[2],controls[4]);

    else
        z_2     = ((1.0+pa.σ)/(ss.λ*pa.β)*(A_cons+(1.0-pa.α)/pa.α*ss.ω*n_full_info))^(1.0/(1.0+pa.σ))
        #Keep the lower number:
        z_upbar = min(z_upbar,z_2);
        z_int_lwbar = max(0.0, (A_cons + pa.ς*ss.ω*(1.0-pa.α)/pa.α)/(ss.λ*pa.β)*(1.0+pa.σ) )^(1.0/(1.0+pa.σ)) # min z for interior n
        weirdinterior = z_int_lwbar<pre_tax_profits_at_nmin/ss.λ && zfoc_at_nmin(z_int_lwbar)>=0 # Boolean for when A is not that big
        if nn_z0 >= pa.ς || weirdinterior

            verbose && println("Case 1: Interior n")
            # 1.1 Evaluating the interior solution:
            # Using bisection method:
            if zfoc_at_nopt(z_int_lwbar)*zfoc_at_nopt(z_upbar)<0
                controls[4] = find_zero(zfoc_at_nopt, (z_int_lwbar,z_upbar), Bisection());
            else
                println("e: ", ss.e," z_2: ", z_2)
                println("z_int_lwbar: ",z_int_lwbar, " z_upbar: ", z_upbar, " zfoc(z_upbar): ", zfoc_at_nopt(z_upbar))
                error("Bisection: signs equal --> Cannot solve.")
            end

            controls[2] = n_opt(controls[4]);
            den_l(controls[2],controls[4]) <= 0.0 && error("Denominator of L is negative or zero.") #Stop when denominator is negative.
            controls[1] = ll_opt(controls[2],controls[4]);
            controls[3] = (pa.χ*controls[1]^(1.0+pa.ψ))/den_p(controls[2],controls[4]);

            max_hamiltonian = objective(controls[1],controls[2],controls[3],controls[4])

            # 1.1.2 Evaluating the corner solution:
            corner_zz = min(z_int_lwbar, pre_tax_profits_at_nmin/ss.λ, z_upbar)
            corner_nn = max(nn_z0,pa.ς)
            den_l(corner_nn,corner_zz) <= 0.0 && error("Denominator of L (corner) is negative or zero.") #Stop when denominator is negative.
            corner_ll = ll_opt(corner_nn,corner_zz)
            corner_pp = (pa.χ*corner_ll^(1.0+pa.ψ))/den_p(corner_nn,corner_zz);

            corner_hamiltonian = objective(corner_ll,corner_nn,corner_pp,corner_zz)

            #3. Defining the solution we are keeping:
            if corner_hamiltonian>max_hamiltonian
            	verbose && println("Case 1.b Corner n and z.")
	                controls[4]=corner_zz
	                controls[2]=corner_nn
	                controls[1]=corner_ll
	                controls[3]=corner_pp
	                max_hamiltonian=corner_hamiltonian
			end
        else
	        verbose && print("Case 2: n at lower bound. ")
	        # 1.2 Evaluate n=nmin
	        z_upbar = min(z_int_lwbar, pre_tax_profits_at_nmin/ss.λ, z_upbar) # z_upbar just for the case \bar m > n_full_info
	        controls[2] = pa.ς
	        # verbose && println("Lower=",zfoc_at_nmin(z_lwbar), "  Upper=", zfoc_at_nmin(z_upbar), "  n=", controls[2])
	        # 1.2.1 Interior z
	        if zfoc_at_nmin(z_lwbar)*zfoc_at_nmin(z_upbar)<0
	            controls[4] = find_zero(zfoc_at_nmin, (z_lwbar,z_upbar), Bisection())
	            if den_l(controls[2],controls[4]) > 0.0
	            	# verbose && print("Interior z feasible. ")
	                controls[1] = ll_opt(controls[2],controls[4]);
	                controls[3] = (pa.χ*controls[1]^(1.0+pa.ψ))/den_p(controls[2],controls[4]);
	                max_hamiltonian = objective(controls[1], controls[2], controls[3], controls[4])
	            end
	        end
	        # 1.2.2 True corner
	        if z_upbar>0.0 && den_l(controls[2], z_upbar)>0.0
	            corner_ll=ll_opt(controls[2], z_upbar);
	            corner_pp = (pa.χ*corner_ll^(1.0+pa.ψ))/den_p(controls[2], z_upbar);
	            corner_hamiltonian = objective(corner_ll, controls[2], corner_pp, z_upbar)
	            if corner_hamiltonian>max_hamiltonian
	            	verbose && println("Case 2.b Corner z.")
	                controls[4]=z_upbar
	                # controls[2]=controls[2]
	                controls[1]=corner_ll
	                controls[3]=corner_pp
	                max_hamiltonian=corner_hamiltonian
	            else
	            	verbose && println("Case 2.a Interior z.")
	            end
	        end

	        # 1.2.3. Evaluating the z-corner solutions (this shouldn't be the answer)
	        corner_zz = 0.0
	        # corner_nn = controls[2];
	        den_l(controls[2], corner_zz) <= 0.0 && error("Denominator of L (corner) is negative or zero.") #Stop when denominator is negative.
	        corner_ll = ll_opt(controls[2], corner_zz)
	        corner_pp = (pa.χ*corner_ll^(1.0+pa.ψ))/den_p(controls[2], corner_zz);
	        corner_hamiltonian = objective(corner_ll, controls[2], corner_pp, corner_zz)
	        if corner_hamiltonian > max_hamiltonian
	            controls[4] = corner_zz
	            # controls[2] = corner_nn;
	            controls[1] = corner_ll
	            controls[3] = corner_pp
	            max_hamiltonian=corner_hamiltonian
	            @warn("Case 2.c Corner n z=0.")
	        end
        end
    end

    #Final Output:
    verbose && println("z = ", controls[4], ",  n = ", controls[2], ",  l = ", controls[1], ",  p = ", controls[3])
    
end

function new_find_controls(θ::Float64, ss::State, pa, verbosebool::Bool=false)
    controls=Array{Float64}(undef,4)
    new_find_controls!(controls, θ::Float64, ss::State, pa, verbosebool)
    return Tuple(controls)
end

function recover_controls!(ctrlvec::Array{Float64}, θvec::Array{Float64}, solvec::Array{Float64})
    (Nspan,~)=size(solvec)

    for j=Nspan:-1:1
      θ     = θvec[j];
      uw    = solvec[j,1];
      μ     = solvec[j,2];
      e     = solvec[j,3];
      ϕ_e   = solvec[j,4];
      λ     = solvec[j,6];
      ω     = solvec[j,8];
      ss    = State(e, uw, ϕ_e, μ, λ, ω);

      
      (ll, nn, pp, zz) = new_find_controls( θ, ss, pa)
      # println("θwControls = ", θ)
      # println("z = ", zz, "n = ", nn, "l = ", ll, "p = ", pp)

      ctrlvec[j,1] = zz;
      ctrlvec[j,2] = nn;
      ctrlvec[j,3] = ll;
      ctrlvec[j,4] = pp;
    end
    nothing
end
