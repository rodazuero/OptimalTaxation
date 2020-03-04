
function find_states!(du,u,pa,θ,ini, par, cons)
    μ     = u[1];
    e     = u[2];

    #Construct state object
    #ss = State(μ,e);
    #ss = State(ini[1], ini[2]);
    #ss = State(e, uw, ϕ_e, μ, λ, ω);
    ss = State(ini[2], 0.0, 0.0, ini[1], 0.0, 0.0);

    #Find optimal controls
    (p) = new_find_controls( ini[3], ss, pa, par, cons);

    h_e   = pa.he( θ, e);
    h_w   = pa.hw( θ, e);
    h_tot = pa.hh( θ, e, p);

    λ     = 1.0;
    du[1] = λ*h_tot;
    du[2] = p;

    nothing
end

function find_states_RKPack!(du,u,pa,θ,ini,par,cons)
    μ     = u[1];
    e     = u[2];

    #Construct state object
    #ss = State(μ,e);
    #ss = State(ini[1], ini[2]);
    #ss = State(e, uw, ϕ_e, μ, λ, ω);
    ss = State(ini[2], 0.0, 0.0, ini[1], 0.0, 0.0);

    #Find optimal controls
    (p) = new_find_controls( ini[3], ss, pa, par, cons);

    h_e   = pa.he( θ, e);
    h_w   = pa.hw( θ, e);
    h_tot = pa.hh( θ, e, p);

    λ     = 1.0;
    du[1] = λ*h_tot;
    du[2] = p;

    nothing
end
