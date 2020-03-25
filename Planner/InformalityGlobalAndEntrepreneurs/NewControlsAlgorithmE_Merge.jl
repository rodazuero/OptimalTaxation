function new_find_controlse(θ::Float64, θe::Float64, sse, pa)
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    #Recover ditributions
    h_e= pa.he(θ, θe);

    #Defining bounds we use in various cases (limits of z):
    #println("θe = ", θe, "μe = ", sse.μe)
    n_full_info  = ((sse.λe*pa.α*θe)/sse.ωfe)^(1.0/(1.0-pa.α));
    z_max   = (1.0/pa.β)^(1.0/pa.σ); #Max possible evasion.
    z_lwbar = eps(); #The smallest possible z.
    z_1     = z_max;
    z_2     = -pa.σ*sse.μe*n_full_info^pa.α/(sse.λe*h_e);
    #Keep the lower number:
    z_upbar = min(z_1,z_2);
    #println("z_lwbar = ", z_lwbar, "z_upbar = ", z_upbar)

    #Defining the functions we are using:
    #fun_ni: entrepreneurs FOC wrt ni
    fun_ni(n)   = ((pa.α*θe*n^(pa.α-1.0)-sse.wie)/pa.δ)^(1.0/pa.γ);
    #n_opt: solving n from planner's Entrepreneurs FOC wrt z
    n_opt(z)    = (-sse.λe*z*h_e/(pa.σ*sse.μe))^(1.0/pa.α);
    #n_opt: Planner's Entrepreneurs FOC wrt n
    fun_nz(z,n) = (sse.λe*pa.α*θe*n^(pa.α-1.0)*h_e - sse.ωfe*h_e + pa.α*sse.μe*n^(pa.α-1.0)*(1.0-pa.β*z^pa.σ) +
                   fun_ni(n)^(1.0-pa.γ)/n^(2.0-pa.α)*pa.α*(1.0-pa.α)*θe/(pa.δ*pa.γ)*(sse.λe*pa.δ*fun_ni(n)^pa.γ+sse.ωie-sse.ωfe)*h_e);
    fun_z(z)    = fun_nz(z,n_opt(z));

    #Hamiltonian:
    objective(ze,ne,nie) = pa.indicator*sse.ue^pa.ϕ*h_e + sse.λe*( θe*ne^pa.α - pa.β*ze^(1.0+pa.σ)/(1.0+pa.σ) - pa.δ/(1.0+pa.γ)*nie^(1.0+pa.γ) - sse.ue)*h_e -
                           sse.ωfe*(ne-nie)*h_e + sse.μe*ne^pa.α*(1.0-pa.β*ze^pa.σ) - sse.ωie*nie*h_e;
       #1.Grid
       Nz = 100; #Number of Zs
       zstep = (z_upbar - z_lwbar)/(Nz-1);
       zgrid = range(z_lwbar, stop = z_upbar, length = Nz);
       #println("Number of z = ", Nz)

       iteration=NaN
       grid_z = Array{Float64}(undef,Nz,4);
       fill!(grid_z,NaN);

       for j = 1:Nz

           #println(j)
           potential_z   = collect(zgrid)[j];
           potential_n   = n_opt(potential_z);
           potential_cpo = fun_z(potential_z);

           grid_z[j,1] = potential_z
           grid_z[j,2] = potential_n
           grid_z[j,3] = potential_cpo
           #To see if there is a change of sign in the derivative:
           if j>1
               grid_z[j,4] = grid_z[j,3]*grid_z[j-1,3]
           else
               grid_z[j,4] = 0
           end
       end

       interior  = grid_z[:,4].>=0
       iteration = findall(x->x==false, interior)

       if all(x->x==true, interior)
           #iteration  = Nz
           #potential_z = grid_z[iteration,1]
           potential_z  = 0.0
           potential_n  = 0.0
           potential_ni = 0.0
           interior_hamiltonian = objective(potential_z,potential_n,potential_ni);
       else
           iteration = iteration[1]
           #println(iteration[1])
           new_lwbar = grid_z[iteration-1,1]
           new_upbar = grid_z[iteration,1]

           if fun_z(new_lwbar)*fun_z(new_upbar)<0
               potential_z = find_zero(fun_z, (new_lwbar,new_upbar), Bisection());
           else
               error("Bisection: signs equal --> Cannot solve.")
           end

          potential_n  = n_opt(potential_z);
          potential_ni = fun_ni(potential_n);

          interior_hamiltonian = objective(potential_z,potential_n,potential_ni);
       end

       #2. The corner solution is:
       corner_zz  = 0.0;
       corner_nn  = 0.0;
       corner_nni = 0.0;

       corner_hamiltonian = objective(corner_zz,corner_nn,corner_nni);

       if corner_hamiltonian > interior_hamiltonian
           zze  = corner_zz;
           nne  = corner_nn;
           nnie = corner_nni;
       else
           zze  = potential_z;
           nne  = potential_n;
           nnie = potential_ni;
       end

       #Final Output:
       println("zze = ", zze, "nne = ", nne, "nnie = ", nnie)
       zze, nne, nnie
   end

   function recover_controlse!(ctrlvec::Array{Float64,2}, θw::Float64 ,θvec::Array{Float64}, solvec::Array{Float64,2})
       (Nspan,~)=size(solvec)

       for j=Nspan:-1:1
           θe      = θvec[j];
           ue      = solvec[j,1];
           μe      = solvec[j,2];
           ye_agg  = solvec[j,3];
           λe      = solvec[j,4];
           le_agg  = solvec[j,5];
           ωfe     = solvec[j,6];
           lie_agg = solvec[j,7];
           ωie     = solvec[j,8];
           wie     = solvec[j,9];
           ϕie     = solvec[j,10];

           sse = StateE(ue, μe, ye_agg, λe, le_agg, ωfe, lie_agg, ωie, wie, ϕie);

           (zze, nne, nnie) = new_find_controlse( θw, θe, sse, pa)

           ctrlvec[j,1] = zze;
           ctrlvec[j,2] = nne;
           ctrlvec[j,3] = nnie;
       end

       nothing
   end
