function new_find_controlse(θ::Float64, θe::Float64, sse, pa)
    #INPUT: states and parameters
    #OUTPUT: optimal controls

    #Recover ditributions
    h_e= pa.he(θ, θe);

    #Defining bounds we use in various cases (limits of z):
    println("θe = ", θe, "μe = ", sse.μe)
    n_full_info  = ((sse.λe*pa.α*θe)/sse.ωfe)^(1.0/(1.0-pa.α));
    z_max   = (1.0/pa.β)^(1.0/pa.σ); #Max possible evasion.
    z_lwbar = eps(); #The smallest possible z.
    z_1     = z_max;
    z_2     = -pa.σ*sse.μe*n_full_info^pa.α/(sse.λe*h_e);
    #Keep the lower number:
    z_upbar = min(z_1,z_2);

    #Defining the functions we are using:
    #fun_ni: entrepreneurs FOC wrt ni
    fun_ni(n)   = ((pa.α*θe*n^(pa.α-1.0)-sse.wie)/pa.δ)^(1.0/pa.γ);
    #n_opt: solving n from planner's Entrepreneurs FOC wrt z
    n_opt(z)    = (-sse.λe*z*h_e/(pa.σ*sse.μe))^(1.0/pa.α);
    #n_opt: Planner's Entrepreneurs FOC wrt n
    fun_nz(z,n) = (sse.λe*pa.α*θe*n^(pa.α-1.0)*h_e - sse.ωfe*h_e + pa.α*sse.μe*n^(pa.α-1.0)*(1.0-pa.β*z^pa.σ)+
                   fun_ni(n)^(1.0-pa.γ)/n^(2.0-pa.α)*pa.α*(1.0-pa.α)*θe/(pa.δ*pa.γ)*(sse.λe*pa.δ*fun_ni(n)^pa.γ+sse.ωie-sse.ωfe)*h_e);
    fun_z(z)    = fun_nz(z,n_opt(z));

    #Hamiltonian:
    objective(ze,ne,nie) = pa.indicator*sse.ue^pa.ϕ*h_e + sse.λe*( θe*ne^pa.α - pa.β*ze^(1.0+pa.σ)/(1.0+pa.σ) - pa.δ/(1.0+pa.γ)*nie^(1.0+pa.γ) - sse.ue)*h_e -
                           sse.ωfe*(ne-nie)*h_e + sse.μe*ne^pa.α*(1.0-pa.β*ze^pa.σ) - sse.ωie*nie*h_e;

    #Defining values for the grid:
    #Nz = 100000; #Number of Zs
    Nz = 10000; #Number of Zs
    zstep = (z_upbar - z_lwbar)/(Nz-1);
    zgrid = range(z_lwbar, stop = z_upbar, length = Nz);
    #println("Number of z = ", Nz)

    #Define a variable counting the number of gridpoints where no solution if found
    global nosol = 0;
    global current_max  = NaN;
    global nne = NaN;
    global zze = NaN;
    global nnie = NaN;

    #1. The interior solution is:
    for j = 1:Nz

        #Solution for f_n:
        potential_z  = zgrid[j];
        potential_n  = n_opt(potential_z);
        potential_ni = fun_ni(potential_n);

        potential_value = objective(potential_z,potential_n,potential_ni);

        #Evaluating the value of the Hamiltonian:
        if isnan(potential_value)
            global nosol = nosol + 1;
        elseif isnan(current_max)
            global current_max = potential_value;
            global nne         = potential_n;
            global zze         = potential_z;
            global nnie        = potential_ni;
        elseif potential_value > current_max
            global current_max = potential_value;
            global nne         = potential_n;
            global zze         = potential_z;
            global nnie        = potential_ni;
        end
    end

    #2. The corner solution is:
    corner_zz  = 0.0;
    corner_nn  = 0.0;
    corner_nni = 0.0;

    corner_hamiltonian = objective(corner_zz,corner_nn,corner_nni);

    if corner_hamiltonian > current_max
        zz  = corner_zz;
        nn  = corner_nn;
        nni = corner_nni;
    else
        zz  = zze;
        nn  = nne;
        nni = nnie;
    end
    #Final Output:
    #println("zze = ", zze, "nne = ", nne, "nnie = ", nnie)
    zz, nn, nni
end

function recover_controlse!(ctrlvec::Array{Float64,2}, θw::Float64 ,θvec::Array{Float64}, solvec::Array{Float64,2})
    (Nspan,~)=size(solvec)

    for j=Nspan:-1:1
        println("j = ", j)
        θe      = θvec[j];
        ue      = solvec[j,1];
        μe      = solvec[j,2];
        ye_agg  = solvec[j,3];
        λe      = solvec[j,4];
        lfe_agg = solvec[j,5];
        ωfe     = solvec[j,6];
        lie_agg = solvec[j,7];
        ωie     = solvec[j,8];
        wie     = solvec[j,9];
        ϕwe     = solvec[j,10];

        sse = StateE(ue, μe, ye_agg, λe, lfe_agg, ωfe, lie_agg, ωie, wie, ϕwe);

        (zze, nne, nnie) = new_find_controlse( θw, θe, sse, pa)

        ctrlvec[j,1] = zze;
        ctrlvec[j,2] = nne;
        ctrlvec[j,3] = nnie;
    end

    nothing
end
