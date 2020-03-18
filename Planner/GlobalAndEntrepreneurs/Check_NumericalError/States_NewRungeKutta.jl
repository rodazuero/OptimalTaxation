
function find_states!(du,u,pa,θ,ini, par, cons)
    μ     = u[1];
    e     = u[2];

    #Find optimal controls
    p = (1.0+par)*cons*(θ-pa.θ_w_lb)^par

    h_e   = pa.he( θ, e);
    h_w   = pa.hw( θ, e);
    h_tot = pa.hh( θ, e, p);

    λ     = 1.0;
    du[1] = λ*h_tot;
    du[2] = p;

    nothing
end
