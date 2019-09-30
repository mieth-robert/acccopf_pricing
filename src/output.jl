function jump_to_df(m, settings, buses, generators, lines, farms; solvetime=0)

    n_buses = length(buses)
    n_lines = length(lines)
    n_generators = length(generators)
    n_farms = length(farms)
    θᵘ = deg2rad(float(theta_u))

    haswind = falses(n_buses)
    for k in 1:n_farms
        haswind[farms[k].bus] = true
    end

    PV = findall(b -> b.kind==:PV, buses)
    PQ = findall(b -> b.kind==:PQ, buses)
    REF = findall(b -> b.kind==:Ref, buses)

    num_approx = 1e-7

    res_buses = collect(1:n_buses)
    res_bustype =  [i in PQ ? "PQ" : (i in PV ? "PV" : "REF") for i in 1:n_buses]

    gp_cc = value.(m[:gp_cc])
    gp_at_bus = [length(buses[i].genids) > 0 ? sum(gp_cc[k] for k in buses[i].genids) : 0 for i in 1:n_buses]
    
    if (settings["run_type"]=="full_cc") || (settings["run_type"]=="gen_cc")
        delta_p_p = shadow_price.(m[:δp])
        delta_p_p_atbus = [length(buses[i].genids) > 0 ? sum(delta_p_p[k] for k in buses[i].genids) : 0 for i in 1:n_buses]
        delta_p_m = shadow_price.(m[:δm])
        delta_p_m_atbus = [length(buses[i].genids) > 0 ? sum(delta_p_m[k] for k in buses[i].genids) : 0 for i in 1:n_buses]
    else
        delta_p_p_atbus = fill(0, n_buses)
        delta_p_m_atbus = fill(0, n_buses)
    end

    delta_p_det_p = shadow_price.(m_cc[:δ_p_det_plus])
    delta_p_det_p_atbus = [length(buses[i].genids) > 0 ? sum(delta_p_det_p[k] for k in buses[i].genids) : 0 for i in 1:n_buses]

    delta_p_det_m = shadow_price.(m_cc[:δ_p_det_min])
    delta_p_det_m_atbus = [length(buses[i].genids) > 0 ? sum(delta_p_det_m[k] for k in buses[i].genids) : 0 for i in 1:n_buses]

    gq_cc = value.(m[:gq_cc])
    gq_at_bus = [length(buses[i].genids) > 0 ? sum(gq_cc[k] for k in buses[i].genids) : 0 for i in 1:n_buses]

    v = value.(m[:v_cc])
    mu_v_p = shadow_price.(m_cc[:μp])
    mu_v_m = shadow_price.(m_cc[:μm])

    if (settings["alpha_mod"] == "split") || (settings["run_type"] == "det")
        alpha_bus = zeros(n_buses)
    else
        alpha_bus = value.(m[:α_by_bus])
    end

    std_v_all = zeros(n_buses)
    if settings["run_type"] == "full_cc"
        std_v_PQ = value.(m[:t_v])
        for (i,b) in enumerate(PQ)
            std_v_all[b] = std_v_PQ[i]
        end
    end
        
    std_q_all = zeros(n_buses)
    if settings["run_type"] == "full_cc"
        std_q_PV = value.(m[:t_qPV])
        for (i,b) in enumerate(PV)
            std_q_all[b] = std_q_PV[i]
        end
    end

    if settings["alpha_mod"] == "split"
        std_p = value.(m[:t_p])
    else
        std_p = value.(m[:α])
    end
    std_p_all = [length(buses[i].genids) > 0 ? sum(std_p[k] for k in buses[i].genids) : 0 for i in 1:n_buses]

    λ_p = shadow_price.(m[:λ_p])
    λ_q = shadow_price.(m[:λ_q])

    γ_all = zeros(n_buses)
    if settings["run_type"] != "det"
        if settings["alpha_mod"] == "split"
            γ = shadow_price.(m[:γ])
            for (i,f) in enumerate(farms)
                γ_all[f.bus] = γ[i]
            end
        else
            γ_all[1] = shadow_price(m[:γ])    
        end
    end


    bus_results = DataFrame(
        bus = res_buses,
        btype = res_bustype,
        haswind = haswind,
        gp = gp_at_bus,
        gq = gq_at_bus,
        v = v,
        alpha = alpha_bus,
        std_v = std_v_all,
        std_q_G = std_q_all,
        std_p_G = std_p_all,
        lambda_p = λ_p,
        lambda_q = λ_q,
        gamma = γ_all, 
        delta_p_plus = delta_p_p_atbus,
        delta_p_minus = delta_p_m_atbus,
        delta_p_det_plus = delta_p_det_p_atbus,
        delta_p_det_minus = delta_p_det_m_atbus,
        mu_v_p = mu_v_p,
        mu_v_m = mu_v_m,
    )


    f_p_ij = value.(m[:pij_cc])
    f_p_ji = value.(m[:pji_cc])
    f_q_ij = value.(m[:qij_cc])
    f_q_ji = value.(m[:qji_cc])

    if settings["run_type"] == "full_cc"
        std_pij = value.(m[:t_f_pij])
        std_pji = value.(m[:t_f_pji])
        std_qij = value.(m[:t_f_qij])
        std_qji = value.(m[:t_f_qji])
        aux_pij = value.(m[:aux_pij])
        aux_pji = value.(m[:aux_pji])
        aux_qij = value.(m[:aux_qij])
        aux_qji = value.(m[:aux_qji])
    else
        std_pij = zeros(n_lines)
        std_pji = zeros(n_lines)
        std_qij = zeros(n_lines)
        std_qji = zeros(n_lines)
        aux_pij = zeros(n_lines)
        aux_pji = zeros(n_lines)
        aux_qij = zeros(n_lines)
        aux_qji = zeros(n_lines)
    end
        
    eta_ij = sum(shadow_price.(m[:η_ij])[:,b] for b in 1:12)
    eta_ji = sum(shadow_price.(m[:η_ji])[:,b] for b in 1:12)

    flow_results = DataFrame(
        f_p_ij = f_p_ij,
        f_p_ji = f_p_ji,
        f_q_ij = f_q_ij,
        f_q_ji = f_q_ji,
        std_pij = std_pij,
        std_pji = std_pji,
        std_qij = std_qij,
        std_qji = std_qji,
        aux_pij = aux_pij,
        aux_pji = aux_pji,
        aux_qij = aux_qij,
        aux_qji = aux_qji,
        eta_ij = eta_ij,
        eta_ji = eta_ji,
    )
          

    objective = [objective_value(m)]
    gen_cost = [value(m[:quad_cost]) + value(m[:linear_cost])]


    if settings["run_type"] == "det"
        var_pen = 0
    else
        var_pen = [value(m[:var_penalty])]
    end

    sum_p_var = sum(value.(m[:t_p]).^2)
    sum_q_var = sum(value.(m[:t_qPV]).^2)
    sum_v_var = sum(value.(m[:t_v]).^2)
    sum_fp_var = sum(value.(m[:t_f_pij]).^2) + sum(value.(m[:t_f_pji]).^2)
    sum_fq_var = sum(value.(m[:t_f_qij]).^2) + sum(value.(m[:t_f_qji]).^2)

    status = termination_status(m)

    objective_results = DataFrame(
        status = status,
        solvetime = solvetime,
        objective = objective,
        gen_cost = gen_cost,
        var_pen = var_pen,
        sum_p_var = sum_p_var,
        sum_q_var = sum_q_var,
        sum_v_var = sum_v_var,
        sum_fp_var = sum_fp_var,
        sum_fq_var = sum_fq_var,
    )


   
    gen = 1:n_generators    
    at_bus = [g.busidx for g in generators]
    alphas = value.(m_cc[:α])
    if settings["alpha_mod"] == "split"
        alpha_results = DataFrame(alphas) 
        farmsbuses = [f.bus for f in farms]
        names!(alpha_results, map(Symbol, farmsbuses))
    else
        alpha_results = DataFrame(alpha = alphas) 
    end
    insertcols!(alpha_results, 1, gen = gen)
    insertcols!(alpha_results, 2, at_bus = at_bus)
    

    return [bus_results, flow_results, objective_results, alpha_results]

end


function save_results(res_dfs, settings, ex_name)
    timestamp = Dates.format(Dates.now(), "yymmdd_HHMM")
    resultdf_names = ["bus_results", "line_results", "objective_results", "alpha_results"]
    savepath = "results/$(ex_name)_$(timestamp)"
    for (name, res_df) in res_dfs
        mkpath("$(savepath)/$(name)")
        for (i,res) in enumerate(resultdf_names)
            CSV.write("$(savepath)/$(name)/$(res).csv", res_df[i])
        end
        settings_str = JSON.json(settings[name])
        open("$(savepath)/$(name)/settings.json", "w") do f
            write(f, settings_str)
        end
    end    
end