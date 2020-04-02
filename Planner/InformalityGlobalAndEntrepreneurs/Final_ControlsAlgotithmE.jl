function new_find_controlse(θe::Float64, sse, pa)
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    #Recover ditributions
    h_e = pa.he(pa.θ_w_ub, θe);

    #Defining definitions we are using:
    n_full_info = ((sse.λe*pa.α*θe)/sse.ωfe)^(1.0/(1.0-pa.α));
    z_max   = (1.0/pa.β)^(1.0/pa.σ); #Max possible evasion.
    z_lwbar = eps(); #The smallest possible z.
    z_1     = z_max;
    z_2     = -pa.σ*sse.μe*n_full_info^pa.α/(sse.λe*h_e);
    z_upbar = min(z_1,z_2); #Keep the lower number as the upper bound
    #println("z_lwbar = ", z_lwbar, "z_upbar = ", z_upbar)

    #Defining the functions we are using:
    #fun_ni: entrepreneurs FOC wrt ni
    fun_ni(n)   = ((pa.α*θe*n^(pa.α-1.0)-sse.wie)/pa.δ)^(1.0/pa.γ);
    #n_opt: solving n from planner's Entrepreneurs FOC wrt z
    n_opt(z)    = (-sse.λe*z*h_e/(pa.σ*sse.μe))^(1.0/pa.α);
    #fun_nz: Planner's Entrepreneurs FOC wrt n
    function fun_nz(z,n)
        if pa.γ == pa.ρ
            pa.α*sse.μe*n^(pa.α-1.0)*(1.0-pa.β*z^pa.σ) + sse.λe*pa.α*θe*n^(pa.α-1.0)*h_e*(1.0 + (1.0-pa.α)/pa.γ*fun_ni(n)/n) - sse.ωfe*h_e;
        else
            (sse.λe*pa.α*θe*n^(pa.α-1.0)*h_e - sse.ωfe*h_e + pa.α*sse.μe*n^(pa.α-1.0)*(1.0-pa.β*z^pa.σ) +
             fun_ni(n)^(1.0-pa.γ)/n^(2.0-pa.α)*pa.α*(1.0-pa.α)*θe/(pa.δ*pa.γ)*(sse.λe*pa.δ*fun_ni(n)^pa.γ+sse.ωie-sse.ωfe)*h_e);
        end
    end
    #fun_z: Planner's Entrepreneurs FOC wrt n, with the optimal n
    fun_z(z) = fun_nz(z,n_opt(z));

    #Hamiltonian:
    objective(ze,ne,nie) = pa.indicator*sse.ue^pa.ϕ*h_e + sse.λe*( θe*ne^pa.α - pa.β*ze^(1.0+pa.σ)/(1.0+pa.σ) - pa.δ/(1.0+pa.γ)*nie^(1.0+pa.γ) - sse.ue)*h_e -
                           sse.ωfe*(ne-nie)*h_e + sse.μe*ne^pa.α*(1.0-pa.β*ze^pa.σ) - sse.ωie*nie*h_e;

    #In the corner solution all variables are zero (we start from z=0):
    current_z     = 0.0;
    current_n     = 0.0;
    current_ni    = 0.0;
    current_ham   = objective(current_z,current_n,current_ni);
    candidate_z   = NaN;
    candidate_n   = NaN;
    candidate_ni  = NaN;
    candidate_ham = NaN;
    previousCPOpositive = false; #Indicates change from positive to negative in the FOC.
    candidate_cpo = Array{Float64,1}

       #1.Grid
       Nz = 100; #Number of Zs
       zstep = (z_upbar - z_lwbar)/(Nz-1);
       zgrid = range(z_lwbar, stop = z_upbar, length = Nz);
       #println("Number of z = ", Nz)

       for j = 1:Nz

           candidate_cpo = fun_z(zgrid[j]);

           if previousCPOpositive == 1
               if candidate_cpo < 0
                   candidate_z   = find_zero(fun_z, (zgrid[j-1],zgrid[j]), Bisection());
                   candidate_n   = n_opt(candidate_z);
                   candidate_ni  = fun_ni(candidate_n);
                   candidate_ham = objective(candidate_z,candidate_n,candidate_ni);
                   if candidate_ham > current_ham
                       current_z   = candidate_z;
                       current_n   = candidate_n;
                       current_ni  = candidate_ni;
                       current_ham = candidate_ham;
                   end #if
                   previousCPOpositive = false;
               end #if
           else
               candidate_cpo > 0 ? previousCPOpositive = true : previousCPOpositive = false;
           end #if
       end #for

       if candidate_cpo>0 #here the value corresponds to the last value in the grid (z_upbound)
           candidate_z   = z_upbar;
           candidate_n   = n_opt(candidate_z);
           candidate_ni  = fun_ni(candidate_n);
           candidate_ham = objective(candidate_z,candidate_n,candidate_ni);
           if candidate_ham > current_ham
               current_z  = candidate_z;
               current_n  = candidate_n;
               current_ni = candidate_ni;
               current_ham = candidate_ham;
           end #if
       end #if

       return current_z, current_n, current_ni;
   end

   function recover_controlse!(ctrlvec::Array{Float64,2},θvec::Array{Float64}, solvec::Array{Float64,2})
       (Nspan,~)=size(solvec)

       for j=Nspan:-1:1
           println(j)
           θe      = θvec[j];
           ue      = solvec[j,1];
           μe      = solvec[j,2];
           ye_agg  = solvec[j,3];
           λe      = solvec[j,4];
           lfe_agg  = solvec[j,5];
           ωfe     = solvec[j,6];
           lie_agg = solvec[j,7];
           ωie     = solvec[j,8];
           wie     = solvec[j,9];
           ϕwe     = solvec[j,10];

           sse = StateE(ue, μe, ye_agg, λe, lfe_agg, ωfe, lie_agg, ωie, wie, ϕwe);

           (zze, nne, nnie) = new_find_controlse( θe, sse, pa)
           println("θeControls = ", θe)
           println("ze = ", zze, "ne = ", nne, "nie = ", nnie)

           ctrlvec[j,1] = zze;
           ctrlvec[j,2] = nne;
           ctrlvec[j,3] = nnie;
       end

       nothing
   end
