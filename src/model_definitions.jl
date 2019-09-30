
function build_ac_opf(thermalLimitscale, generators, buses, lines, farms, settings)

    ϵ = settings["ϵ"]
    theta_u = settings["theta_u"]

    n_buses = length(buses)
    n_lines = length(lines)
    n_generators = length(generators)
    n_farms = length(farms)
    θᵘ = deg2rad(float(theta_u))

    m = Model(with_optimizer(Ipopt.Optimizer))

    PV = findall(b -> b.kind==:PV, buses)
    PQ = findall(b -> b.kind==:PQ, buses)
    REF = findall(b -> b.kind==:Ref, buses)
    @assert length(REF) == 1

    s = sqrt(sum([f.σ^2 for f in farms]))
    s_sq = s^2
    z_ϵ = quantile(Normal(0, 1), 1-ϵ)
    total_deviation = z_ϵ * s

    println(">>> Total Deviaton = $(total_deviation)")

    # Continuous Variables
    @variable(m, pij[1:n_lines])
    @variable(m, qij[1:n_lines])
    @variable(m, pji[1:n_lines])
    @variable(m, qji[1:n_lines])
    @variable(m, θ[1:n_buses])
    @variable(m, -θᵘ <= θij[1:n_lines] <= θᵘ)
    @variable(m, gp[1:n_generators])
    @variable(m, gq[1:n_generators])
    @variable(m, r[1:n_generators] >= 0)
    @variable(m, v[1:n_buses] >= 0)
    @variable(m, γ[1:n_farms] >= 0, start=0.5) # Only a decision variable at buses with farms

    @constraint(m, sum(r) >= total_deviation)

    #Objective
    @NLobjective(m, Min, sum(generators[i].pi1 * gp[i]^2 + generators[i].pi2 * gp[i] + generators[i].pi3 for i in 1:n_generators))

    #AC Equations
    @NLconstraint(m, [i=1:n_lines], θij[i] == θ[lines[i].head] - θ[lines[i].tail])
    @NLconstraint(m, pij_constr[i=1:n_lines], pij[i] == (1/lines[i].ratio)^2 * lines[i].γ * v[lines[i].head]^2 - (1/lines[i].ratio)*lines[i].γ * v[lines[i].head]*v[lines[i].tail] *cos(θ[lines[i].head] - θ[lines[i].tail]) - (1/lines[i].ratio)*lines[i].β * v[lines[i].head]*v[lines[i].tail] *sin(θ[lines[i].head] - θ[lines[i].tail]))
    @NLconstraint(m, qij_constr[i=1:n_lines], qij[i] == -(1/lines[i].ratio)^2*(lines[i].β + lines[i].b_charge/2)* v[lines[i].head]^2 + (1/lines[i].ratio)*lines[i].β * v[lines[i].head]*v[lines[i].tail] *cos(θ[lines[i].head] - θ[lines[i].tail]) - (1/lines[i].ratio)*lines[i].γ * v[lines[i].head]*v[lines[i].tail] *sin(θ[lines[i].head] - θ[lines[i].tail]))
    # flow leaving j towards i
    @NLconstraint(m, pji_constr[i=1:n_lines], pji[i] ==  lines[i].γ * v[lines[i].tail]^2 - (1/lines[i].ratio) * lines[i].γ * v[lines[i].head]*v[lines[i].tail] *cos(θ[lines[i].head] - θ[lines[i].tail]) + (1/lines[i].ratio) * lines[i].β * v[lines[i].head]*v[lines[i].tail] *sin(θ[lines[i].head] - θ[lines[i].tail]))
    @NLconstraint(m, qji_constr[i=1:n_lines], qji[i] == -(lines[i].β + lines[i].b_charge/2) * v[lines[i].tail]^2 +(1/lines[i].ratio) * lines[i].β * v[lines[i].head]*v[lines[i].tail] *cos(θ[lines[i].head] - θ[lines[i].tail]) + (1/lines[i].ratio) * lines[i].γ * v[lines[i].head]*v[lines[i].tail] * sin(θ[lines[i].head] - θ[lines[i].tail]))

    # Line limits
    for i in 1:n_lines
        if lines[i].u > 0
            @NLconstraint(m, pij[i]^2 + qij[i]^2 <= lines[i].u^2)
            @NLconstraint(m, pji[i]^2 + qji[i]^2 <= lines[i].u^2)
        end
    end

    # Flow balance
    farms_at_bus = []
    for i in 1:n_farms
        farms_at_i = [f for f in farms if f.bus == i]
        push!(farms_at_bus, farms_at_i)
    end
    @expression(m, sum_pgen[i=1:n_buses], sum(gp[k] for k in buses[i].genids) + (buses[i].farmids == [] ? 0 : sum(farms[k].μ for k in buses[i].farmids)))
    @expression(m, sum_qgen[i=1:n_buses], sum(gq[k] for k in buses[i].genids) + (buses[i].farmids == [] ? 0 : sum(γ[k] * farms[k].μ for k in buses[i].farmids)))
    for i in 1:n_buses
        @constraint(m, sum_pgen[i] - buses[i].Pd - buses[i].Gs * v[i]^2 == sum(pij[k] for k in buses[i].outlist) + sum(pji[k] for k in buses[i].inlist))
        @constraint(m, sum_qgen[i] - buses[i].Qd + buses[i].Bs * v[i]^2 == sum(qij[k] for k in buses[i].outlist) + sum(qji[k] for k in buses[i].inlist))
    end
    @constraint(m, θ[REF] .== 0)


    # Operational limits with reserves
    @constraint(m, c_gen1[i=1:n_generators], gp[i] - r[i] >= generators[i].Pgmin)
    @constraint(m, c_gen2[i=1:n_generators], gp[i] + r[i] <= generators[i].Pgmax)
    @constraint(m, c_gen3[i=1:n_generators], gq[i] >= generators[i].Qgmin)
    @constraint(m, c_gen4[i=1:n_generators], gq[i] <= generators[i].Qgmax)

    @constraint(m, c_v1[i=1:n_buses], v[i] <= buses[i].Vmax)
    @constraint(m, c_v2[i=1:n_buses], v[i] >= buses[i].Vmin)

    return m 
end


function build_ac_cc_opf(generators, buses, lines, farms, settings, m_det; print_output=true)

    a1 = [1, 1, 0.2679, -0.2679, -1, -1, -1, -1, -0.2679, 0.2679, 1, 1]
    a2 = [0.2679, 1, 1, 1, 1, 0.2679, -0.2679, -1, -1, -1, -1, -0.2679]
    a3 = -1 .* [-1, -1.366, -1, -1, -1.366, -1, -1, -1.366, -1, -1, -1.366, -1]

    ϵ = settings["ϵ"]
    theta_u = settings["theta_u"]
    if "alpha_mod" in keys(settings)
        alpha_mod = settings["alpha_mod"]
    else
        alpha_mod = "single"
    end
    if "var_penalty" in keys(settings)
        var_penalty = settings["var_penalty"]
    else
        var_penalty = Dict("p_G" => 0, "p_Q" => 0, "v" => 0, "f_p" => 0, "f_q" => 0)
    end
    if "run_type" in keys(settings)
        run_type = settings["run_type"]
    else
        run_type = "full_cc"
    end

    n_buses = length(buses)
    n_lines = length(lines)
    n_generators = length(generators)
    n_farms = length(farms)
    θᵘ = deg2rad(float(theta_u))

    PV = findall(b -> b.kind==:PV, buses)
    PQ = findall(b -> b.kind==:PQ, buses)
    REF = findall(b -> b.kind==:Ref, buses)
    @assert length(REF) == 1

    haswind = falses(n_buses)
    for k in 1:n_farms
        haswind[farms[k].bus] = true
    end

    s = sqrt(sum([f.σ^2 for f in farms]))
    s_sq = s^2
    z_ϵ = quantile(Normal(0, 1), 1-ϵ)
    total_deviation = s * z_ϵ
    Σ = diagm(0 => [f.σ^2 for f in farms])
    Σ_rt = Σ^(1/2)

    gp_initial = value.(m_det[:gp])
    gq_initial = value.(m_det[:gq])
    v_initial = value.(m_det[:v])
    θ_initial = value.(m_det[:θ])
    γ_initial = value.(m_det[:γ])

    # linearized map from (v,θ) to flows
    x = [m_det[:v]; m_det[:θ]]
    flows = [m_det[:pij];m_det[:pji];m_det[:qij];m_det[:qji]]
    flows_constr = [m_det[:pij_constr];m_det[:pji_constr];m_det[:qij_constr];m_det[:qji_constr]]

    vθ_to_flows_const, vθ_to_flows_J = extract_jacobian(m_det, flows_constr, x)
    vθ_to_flows_J *= -1 # JuMP moves lhs to rhs, model is defined with constants on lhs
    vθ_to_flows_const -= value.(flows)
    vθ_to_flows_const *= -1
    vθ_to_flows_const -= vθ_to_flows_J * value.(x)

    # flows ≈ vθ_to_flows_const + vθ_to_flows_J*(x)

    # now set up map from (v,θ) to nodal (p,q)
    vθ_to_node_p_J = zeros(n_buses, 2 * n_buses)
    vθ_to_node_p_const = zeros(n_buses)
    vθ_to_node_q_J = zeros(n_buses, 2 * n_buses)
    vθ_to_node_q_const = zeros(n_buses)
    for i=1:n_buses
        # p = Pd + Gs*v[i]^2 + sum{pij[k], k in buses[i].outlist} + sum{pji[k], k in buses[i].inlist}
        vθ_to_node_p_const[i] = buses[i].Pd
        vθ_to_node_p_const[i] -= buses[i].Gs*v_initial[i]^2
        vθ_to_node_p_J[i,i] += 2*buses[i].Gs*v_initial[i]
        for k in buses[i].outlist
            vθ_to_node_p_J[i, :] += vθ_to_flows_J[k, :] # pij
            vθ_to_node_p_const[i] += vθ_to_flows_const[k]
        end
        for k in buses[i].inlist
            vθ_to_node_p_J[i,:] += vθ_to_flows_J[k + n_lines, :] # pji
            vθ_to_node_p_const[i] += vθ_to_flows_const[k + n_lines]
        end

        vθ_to_node_q_const[i] = buses[i].Qd
        vθ_to_node_q_const[i] += buses[i].Bs*v_initial[i]^2
        vθ_to_node_q_J[i,i] -= 2*buses[i].Bs*v_initial[i]
        for k in buses[i].outlist
            vθ_to_node_q_J[i,:] += vθ_to_flows_J[k + 2 * n_lines, :] # qij
            vθ_to_node_q_const[i] += vθ_to_flows_const[k + 2 * n_lines]
        end
        for k in buses[i].inlist
            vθ_to_node_q_J[i,:] += vθ_to_flows_J[k + 3 * n_lines, :] # qji
            vθ_to_node_q_const[i] += vθ_to_flows_const[k + 3 * n_lines]
        end
    end

    # p ≈ vθ_to_node_p_const + vθ_to_node_p_J*(x)
    # q ≈ vθ_to_node_q_const + vθ_to_node_q_J*(x)

    # Controllable Generation
    pgen = [sum(Float64[gp_initial[k] for k in buses[i].genids]) for i in 1:n_buses]
    qgen = [sum(Float64[gq_initial[k] for k in buses[i].genids]) for i in 1:n_buses]

    # Total generation based on controlled PF
    for k in 1:n_farms
        pgen[farms[k].bus] += farms[k].μ
        qgen[farms[k].bus] += γ_initial[k]*farms[k].μ
    end

    vθ_to_node_pq_J = [ vθ_to_node_p_J
                        vθ_to_node_q_J ]
    vθ_to_node_pq_const = [ vθ_to_node_p_const
                            vθ_to_node_q_const ]


    #=
    [ y ] = [ A B ][ η ] + [ c₁ ]
    [ z ]   [ C D ][ β ] + [ c₂ ]
    where η = [ v_pq, θ_pq, θ_pv ]
          β = [ v_pv, v_θv, θ_θv ]
          y = [ p_pq, p_pv, q_pq ]
          z = [ p_θv, q_pv, q_θv ].

    y collects genertor reactions depending on control policy decisions
    z collects all implicit generator reactions

    Want to express η as a function of y, β,
        so η = A^{-1}(y - B*β - c₁)
    Then z = C*η + D*β + c₂ is easy to get as well.
    =#

    y_rows = vcat(PQ, PV, PQ .+ n_buses)
    z_rows = vcat(REF, PV .+ n_buses, REF .+ n_buses)
    η_cols = vcat(PQ, PQ .+ n_buses, PV .+ n_buses)
    β_cols = vcat(PV, REF, REF .+ n_buses)

    A = sparse(vθ_to_node_pq_J[y_rows, η_cols])
    B = sparse(vθ_to_node_pq_J[y_rows, β_cols])
    C = sparse(vθ_to_node_pq_J[z_rows, η_cols])
    D = sparse(vθ_to_node_pq_J[z_rows, β_cols])
    c₁ = vθ_to_node_pq_const[y_rows]
    c₂ = vθ_to_node_pq_const[z_rows]

    Ainv = A\diagm(0 => ones(size(A,1)))
    Ainv[abs.(Ainv) .<= 1e-6] .= 0.0

    Ainvsp = sparse(Ainv)

    # Cost matrix
    Cost_mat = zeros(n_generators, n_generators) # quadratic cost matrix
    for i in 1:n_generators
        Cost_mat[i,i] = generators[i].pi1
    end
    C_rt = Cost_mat^(1/2)

    if print_output
        parlog = 1
    else
        parlog = 0
    end

    m_cc = Model(with_optimizer(Mosek.Optimizer, MSK_IPAR_LOG=parlog))

    @variable(m_cc, θ_cc[1:n_buses])
    @variable(m_cc, gp_cc[1:n_generators])
    @variable(m_cc, gq_cc[1:n_generators])
    @variable(m_cc, v_cc[1:n_buses])
    @variable(m_cc, r_cc[1:n_generators] >= 0)
    @variable(m_cc, pij_cc[1:n_lines])
    @variable(m_cc, qij_cc[1:n_lines])
    @variable(m_cc, pji_cc[1:n_lines])
    @variable(m_cc, qji_cc[1:n_lines])
    @variable(m_cc, γ_cc[1:n_farms] >= 0) 
    if alpha_mod == "split"
        @variable(m_cc, 0 <= α[1:n_generators, 1:n_farms] <= 1)
    else
        @variable(m_cc, α[1:n_generators] >= 0)
    end

    @constraint(m_cc, [i=1:n_generators], generators[i].Qgmin <= gq_cc[i] <= generators[i].Qgmax)

    # Linearization Constraints Deterministic
    @expression(m_cc, sum_pgen[i=1:n_buses], sum(gp_cc[k] for k in buses[i].genids) + (buses[i].farmids == [] ? 0 : sum(farms[k].μ for k in buses[i].farmids)))
    @expression(m_cc, sum_qgen[i=1:n_buses], sum(gq_cc[k] for k in buses[i].genids) + (buses[i].farmids == [] ? 0 : sum(γ_cc[k] * farms[k].μ for k in buses[i].farmids)))
    # p ≈ vθ_to_node_p_const + vθ_to_node_p_J*(x)
    # q ≈ vθ_to_node_q_const + vθ_to_node_q_J*(x)
    x_cc = vcat(v_cc, θ_cc)
    @constraint(m_cc, λ_p, sum_pgen .== vθ_to_node_p_const + vθ_to_node_p_J*x_cc)
    @constraint(m_cc, λ_q, sum_qgen .== vθ_to_node_q_const + vθ_to_node_q_J*x_cc)

    @constraint(m_cc, sum(r_cc) >= total_deviation)

    # Operational line limits with reserves
    @constraint(m_cc, δ_p_det_plus[i=1:n_generators], gp_cc[i] - r_cc[i] >= generators[i].Pgmin)
    @constraint(m_cc, δ_p_det_min[i=1:n_generators], gp_cc[i] + r_cc[i] <= generators[i].Pgmax)

    # [pij,pji,qij,qji] ≈ vθ_to_flows_const + vθ_to_flows_J*(x)
    flows_cc = vcat(pij_cc, pji_cc, qij_cc, qji_cc)
    @constraint(m_cc, flows_cc .== vθ_to_flows_const + vθ_to_flows_J * x_cc)

    # Deterministic line limits
    if (run_type == "det") || (run_type == "gen_cc")
        @constraint(m_cc, η_ij[i=1:n_lines, c=1:12], a1[c]*pij_cc[i] + a2[c]*qij_cc[i] <= a3[c]*lines[i].u)
        @constraint(m_cc, η_ji[i=1:n_lines, c=1:12], a1[c]*pji_cc[i] + a2[c]*qji_cc[i] <= a3[c]*lines[i].u)    
    else
        for i in 1:n_lines
            if lines[i].u > 0
                @constraint(m_cc, [lines[i].u, pij_cc[i], qij_cc[i]] in SecondOrderCone())
                @constraint(m_cc, [lines[i].u, pji_cc[i], qji_cc[i]] in SecondOrderCone())
            end
        end
    end
    

    # Slack bus
    @constraint(m_cc, θ_cc[REF] .== 0)

    ## Probabilistic Response
    if (run_type == "full_cc") || (run_type == "gen_cc")
        if alpha_mod == "split"
            @constraint(m_cc, γ[j=1:n_farms], sum(α[i,j] for i in 1:n_generators) == 1)
            @expression(m_cc, Abus[i=1:n_buses, j=1:n_farms], length(buses[i].genids) > 0 ? sum(α[k,j] for k in buses[i].genids) : 0.0)
        else
            @constraint(m_cc, γ, sum(α[i] for i in 1:n_generators) == 1.)
            @expression(m_cc, α_by_bus[i=1:n_buses], length(buses[i].genids) > 0 ? sum(α[k] for k in buses[i].genids) : 0.0)
        end
    end

    if run_type == "full_cc"
        # Farm-to-node-map
        F = zeros(n_buses, n_farms)
        for i in 1:n_buses
            for j in buses[i].farmids
                F[i,j] = 1
            end
        end
        Fg =  F[:,1] .* γ_cc[1]
        for i in 2:n_farms
            Fg = hcat(Fg, F[:,i] .* γ_cc[i])
        end

        # Voltage sensitivity depending on alpha and gamma

        # first solve for Δη = Δ[v_pq,θ_pq,θ_pv] given Δy = Δ[p_pq, p_pv, q_pq] and Δβ = Δ[v_pv,v_θv,θ_θv].
        # Actually Δβ is zero since they're fixed, don't respond.
        # so Δη = A^{-1}(Δy)

        e_f = ones(n_farms)
        PQPV = vcat(PQ, PV)
        # Only voltage at PQ nodes is senstive to the uncertainty
        VSP = Ainvsp[1:length(PQ), 1:length(PQ)+length(PV)] # Voltage sensitivity to p
        VSQ = Ainvsp[1:length(PQ), length(PQ)+length(PV)+1:end] # Voltage sensitvity to q
        if alpha_mod == "split"
            J_v_α_γ = VSP * (F[PQPV, :] - Abus[PQPV,:]) + VSQ * Fg[PQ, :]
        else
            J_v_α_γ = VSP * (F[PQPV, :] - α_by_bus[PQPV]*e_f') + VSQ * Fg[PQ, :]
        end

        @variable(m_cc, t_v[1:length(PQ)] >= 0) # Standar Deviation (for cc formulation)
        # @variable(m_cc, t_v_sq[1:length(PQ)] >= 0) # Variance (for Varaince penalty)

        t_v_all_buses = []
        for i in 1:n_buses
            if i in PQ
                b = findall(b -> b==i, PQ)[1]
                push!(t_v_all_buses, AffExpr(0.0, t_v[b] => 1.0))
            else
                push!(t_v_all_buses, AffExpr(0.0))
            end
        end
        for i in 1:length(PQ)
            @constraint(m_cc, vec(vcat(t_v[i], (J_v_α_γ[i, :]' * Σ_rt)')) in SecondOrderCone())
            # @constraint(m_cc, vec(vcat(t_v_sq[i], 0.5, (J_v_α_γ[i, :]' * Σ_rt)')) in RotatedSecondOrderCone())
        end        
        # Voltage chance constraints, t is Stdv(v)
        @constraint(m_cc, μp[i=1:n_buses], v_cc[i] + z_ϵ * t_v_all_buses[i] <= buses[i].Vmax)
        @constraint(m_cc, μm[i=1:n_buses], -v_cc[i] + z_ϵ * t_v_all_buses[i] <= -buses[i].Vmin)



        # Now the same for θs at PQ and PV buses
        θSP_PQ = Ainvsp[length(PQ)+1:2*length(PQ), 1:length(PQ)+length(PV)] # theta sensitivity to p
        θSQ_PQ = Ainvsp[length(PQ)+1:2*length(PQ), length(PQ)+length(PV)+1:end] # theta sensitivity to q

        θSP_PV = Ainvsp[2*length(PQ)+1:end, 1:length(PQ)+length(PV)] # theta sensitivity to p
        θSQ_PV = Ainvsp[2*length(PQ)+1:end, length(PQ)+length(PV)+1:end] # theta sensitivity to q
        
        if alpha_mod == "split"
            J_θPQ_α_γ = θSP_PQ * (F[PQPV, :] - Abus[PQPV,:]) + θSQ_PQ * Fg[PQ, :]
            J_θPV_α_γ = θSP_PV * (F[PQPV, :] - Abus[PQPV,:]) + θSQ_PV * Fg[PQ, :]
        else
            J_θPQ_α_γ = θSP_PQ * (F[PQPV, :] - α_by_bus[PQPV]*e_f') + θSQ_PQ * Fg[PQ, :]    
            J_θPV_α_γ = θSP_PV * (F[PQPV, :] - α_by_bus[PQPV]*e_f') + θSQ_PV * Fg[PQ, :]
        end
        

        θSP_PV = Ainvsp[2*length(PQ)+1:end, 1:length(PQ)+length(PV)] # theta sensitivity to p
        θSQ_PV = Ainvsp[2*length(PQ)+1:end, length(PQ)+length(PV)+1:end] # theta sensitivity to q
        

        @variable(m_cc, t_θ[1:n_lines] >= 0)
        for i in 1:n_lines
            ih = lines[i].head
            if ih in PQ
                h = findall(b -> b==ih, PQ)[1]
                J_θ_h = J_θPQ_α_γ[h, :] 
            elseif ih in PV
                h = findall(b -> b==ih, PV)[1]
                J_θ_h = J_θPV_α_γ[h, :] 
            elseif ih in REF
                J_θ_h = zeros(n_farms) # REF bus θ is insensitive to the uncertainty
            end

            it = lines[i].tail
            if it in PQ
                t = findall(b -> b==it, PQ)[1]
                J_θ_t = J_θPQ_α_γ[t, :] 
            elseif it in PV
                t = findall(b -> b==it, PV)[1]
                J_θ_t = J_θPV_α_γ[t, :] 
            elseif it in REF
                J_θ_t = zeros(n_farms) 
            end

            @constraint(m_cc, vec(vcat(t_θ[i], ((J_θ_h - J_θ_t)' * Σ_rt)')) in SecondOrderCone())
        end 
        @constraint(m_cc, νp[i=1:n_lines], (θ_cc[lines[i].head] - θ_cc[lines[i].tail]) + z_ϵ * t_θ[i] <= θᵘ)
        @constraint(m_cc, νm[i=1:n_buses], -(θ_cc[lines[i].head] - θ_cc[lines[i].tail]) + z_ϵ * t_θ[i] <= θᵘ)


        # now compute how the generators need to respond other than α*Ω
        # # Δz = [ p_θv, q_pv, q_θv ]
        # Δz = C*Δη

        # # Active Power at slack 
        # # For now we assume infinite generator at reference 
        # @constraint(m_cc, α_by_bus[REF[1]] == 0)
        # # We must have ω_by_bus[REF] - α_by_bus[REF]Ω + ??? == Δz[1]
        # # So the reference bus adjusts its production to balance the (linearized) system
        # # --> ??? = Δz[1] + α_by_bus[REF]Ω - ω_by_bus[REF]
        # # So the extra production at the reference bus is Δz[1] - ω_by_bus[REF]
        # @variable(m_cc, t_pREF >= 0)
        # J_pREF = C[1,1:length(PQ)]' * J_v_α_γ + C[1,length(PQ)+1:2*length(PQ)]' * J_θPQ_α_γ + C[1, 2*length(PQ)+1:end]' * J_θPV_α_γ
        # J_pREF -= F[REF,:]
        # @constraint(m_cc, vec(vcat(t_pREF, J_pREF * Σ_rt)) in SecondOrderCone())
        # @assert length(buses[REF[1]].genids) == 1 # If this is not true some node-internal error distributuin has to be found
        # REFGEN = buses[REF[1]].genids[1]
        # @constraint(m_cc, gp_cc[REFGEN] + z_ϵ * t_pREF <= generators[REFGEN].Pgmax)
        # @constraint(m_cc, -gp_cc[REFGEN] + z_ϵ * t_pREF <= -generators[REFGEN].Pgmin)

        # Reactive Power at Non-Slack
        @variable(m_cc, t_qPV[1:length(PV)] >= 0)
        J_qPV = C[2:1+length(PV),1:length(PQ)] * J_v_α_γ + C[2:1+length(PV),length(PQ)+1:2*length(PQ)] * J_θPQ_α_γ + C[2:1+length(PV), 2*length(PQ)+1:end]' * J_θPV_α_γ
        @constraint(m_cc, [i=1:length(PV)], vec(vcat(t_qPV[i], (J_qPV[i, :]' * Σ_rt)')) in SecondOrderCone())

        δ_q_cc_plus = []
        δ_q_cc_min = []

        for (i,b) in enumerate(PV)
            @assert length(buses[b].genids) == 1 # see above, PV needs to have at least one generator to be a PV bus
            g = buses[b].genids[1] 
            dqp = @constraint(m_cc, gq_cc[g] + z_ϵ * t_qPV[i] <= generators[g].Qgmax)
            dqm = @constraint(m_cc, -gq_cc[g] + z_ϵ * t_qPV[i] <= -generators[g].Qgmin)
            push!(δ_q_cc_plus, dqp)
            push!(δ_q_cc_min, dqm)
        end

        # # Reactive Power at slack 
        # # For now we assume infinite generator at reference 
        # @variable(m_cc, t_qREF >= 0)
        # J_qREF = C[end,1:length(PQ)]' * J_v_α_γ + C[end,length(PQ)+1:2*length(PQ)]' * J_θPQ_α_γ + C[end, 2*length(PQ)+1:end]' * J_θPV_α_γ
        # @constraint(m_cc, vec(vcat(t_qREF, J_qREF * Σ_rt)) in SecondOrderCone())
        # @constraint(m_cc, gq_cc[REFGEN] + z_ϵ * t_qREF <= generators[REFGEN].Qgmax)
        # @constraint(m_cc, -gq_cc[REFGEN] + z_ϵ * t_qREF <= -generators[REFGEN].Qgmin)


        # Done with generators, now the lines
        # Δ[pij,pji,qij,qji] ≈ vθ_to_flows_J*(Δx)

        J_pij = vθ_to_flows_J[1:n_lines, PQ] * J_v_α_γ              + vθ_to_flows_J[1:n_lines, PQ .+ n_buses] * J_θPQ_α_γ             + vθ_to_flows_J[1:n_lines, PV .+ n_buses] * J_θPV_α_γ
        J_pji = vθ_to_flows_J[n_lines+1:2*n_lines, PQ] * J_v_α_γ    + vθ_to_flows_J[n_lines+1:2*n_lines, PQ .+ n_buses] * J_θPQ_α_γ   + vθ_to_flows_J[n_lines+1:2*n_lines, PV .+ n_buses] * J_θPV_α_γ
        J_qij = vθ_to_flows_J[2*n_lines+1:3*n_lines, PQ] * J_v_α_γ  + vθ_to_flows_J[2*n_lines+1:3*n_lines, PQ .+ n_buses] * J_θPQ_α_γ + vθ_to_flows_J[2*n_lines+1:3*n_lines, PV .+ n_buses] * J_θPV_α_γ
        J_qji = vθ_to_flows_J[3*n_lines+1:end, PQ] * J_v_α_γ        + vθ_to_flows_J[3*n_lines+1:end, PQ .+ n_buses] * J_θPQ_α_γ       + vθ_to_flows_J[3*n_lines+1:end, PV .+ n_buses] * J_θPV_α_γ

        @variable(m_cc, t_f_pij[1:n_lines] >= 0)
        @variable(m_cc, t_f_pji[1:n_lines] >= 0)
        @variable(m_cc, t_f_qij[1:n_lines] >= 0)
        @variable(m_cc, t_f_qji[1:n_lines] >= 0)

        for i in 1:n_lines
            @constraint(m_cc, vec(vcat(t_f_pij[i], (J_pij[i, :]' * Σ_rt)')) in SecondOrderCone())
            @constraint(m_cc, vec(vcat(t_f_pji[i], (J_pji[i, :]' * Σ_rt)')) in SecondOrderCone())
            @constraint(m_cc, vec(vcat(t_f_qij[i], (J_qij[i, :]' * Σ_rt)')) in SecondOrderCone())
            @constraint(m_cc, vec(vcat(t_f_qji[i], (J_qji[i, :]' * Σ_rt)')) in SecondOrderCone())
        end



        # Inner approximation of quadratic line constraints
        @variable(m_cc, aux_pij[1:n_lines] >= 0)
        @variable(m_cc, aux_pji[1:n_lines] >= 0)
        @variable(m_cc, aux_qij[1:n_lines] >= 0)
        @variable(m_cc, aux_qji[1:n_lines] >= 0)

        # apx_fac = 2.5
        apx_fac = 1
        z_inner_left = quantile(Normal(0, 1), ϵ/apx_fac)
        z_inner_right = quantile(Normal(0, 1), 1-(ϵ/apx_fac))
        z_inner_twoside = quantile(Normal(0, 1), 1-(ϵ/(2*apx_fac)))
        
        @constraint(m_cc, [i=1:n_lines], -aux_pij[i] - pij_cc[i] <= z_inner_left * t_f_pij[i])
        @constraint(m_cc, [i=1:n_lines], -aux_pji[i] - pji_cc[i] <= z_inner_left * t_f_pji[i])
        @constraint(m_cc, [i=1:n_lines], -aux_qij[i] - qij_cc[i] <= z_inner_left * t_f_qij[i])
        @constraint(m_cc, [i=1:n_lines], -aux_qji[i] - qji_cc[i] <= z_inner_left * t_f_qji[i])

        @constraint(m_cc, [i=1:n_lines], aux_pij[i] - pij_cc[i] >= z_inner_right * t_f_pij[i])
        @constraint(m_cc, [i=1:n_lines], aux_pji[i] - pji_cc[i] >= z_inner_right * t_f_pji[i])
        @constraint(m_cc, [i=1:n_lines], aux_qij[i] - qij_cc[i] >= z_inner_right * t_f_qij[i])
        @constraint(m_cc, [i=1:n_lines], aux_qji[i] - qji_cc[i] >= z_inner_right * t_f_qji[i])

        @constraint(m_cc, [i=1:n_lines], aux_pij[i] >= z_inner_twoside * t_f_pij[i])
        @constraint(m_cc, [i=1:n_lines], aux_pji[i] >= z_inner_twoside * t_f_pji[i])
        @constraint(m_cc, [i=1:n_lines], aux_qij[i] >= z_inner_twoside * t_f_qij[i])
        @constraint(m_cc, [i=1:n_lines], aux_qji[i] >= z_inner_twoside * t_f_qji[i])

        
        # for i in 1:n_lines
        #     if lines[i].u > 0
        #         @constraint(m_cc, [lines[i].u, aux_pij[i], aux_qij[i]] in SecondOrderCone())
        #         @constraint(m_cc, [lines[i].u, aux_pji[i], aux_qji[i]] in SecondOrderCone())
        #     end
        # end

        # @constraint(m_cc, η_ij[i=1:n_lines], vec(vcat(lines[i].u, aux_pij[i], aux_qij[i])) in SecondOrderCone())
        # @constraint(m_cc, η_ji[i=1:n_lines], vec(vcat(lines[i].u, aux_pji[i], aux_qji[i])) in SecondOrderCone())

        # Polygonal approximation to recover dual variable
        @constraint(m_cc, η_ij[i=1:n_lines, c=1:12], a1[c]*aux_pij[i] + a2[c]*aux_qij[i] <= a3[c]*lines[i].u)
        @constraint(m_cc, η_ji[i=1:n_lines, c=1:12], a1[c]*aux_pji[i] + a2[c]*aux_qji[i] <= a3[c]*lines[i].u)
    else

        # Deterministic Voltage Constraints
        @constraint(m_cc, μp[i=1:n_buses], v_cc[i] <= buses[i].Vmax)
        @constraint(m_cc, μm[i=1:n_buses], -v_cc[i] <= -buses[i].Vmin)


        @expression(m_cc, t_v[i=1:length(PQ)], 0)
        @expression(m_cc, t_qPV[i=1:length(PV)], 0)
        @expression(m_cc, t_f_pji[i=1:n_lines], 0)
        @expression(m_cc, t_f_pij[i=1:n_lines], 0)
        @expression(m_cc, t_f_qji[i=1:n_lines], 0)
        @expression(m_cc, t_f_qij[i=1:n_lines], 0)
    end


    # Active Power at Non-Slack
    # no_slack_gen = setdiff(collect(1:n_generators), findall(g -> g.busidx == REF[1], generators))
    no_slack_gen = collect(1:n_generators)
    if (run_type == "full_cc") || (run_type == "gen_cc")    
        @variable(m_cc, t_p[i=no_slack_gen] >= 0)
        if alpha_mod == "split"
            @constraint(m_cc, [i=no_slack_gen], vec(vcat(t_p[i], (α[i,:]' * Σ_rt)' )) in SecondOrderCone())
        else
            @constraint(m_cc, [i=no_slack_gen], t_p[i] == α[i]*s)
        end
        @constraint(m_cc, δp[i=no_slack_gen], gp_cc[i] + z_ϵ*t_p[i] <= generators[i].Pgmax)
        @constraint(m_cc, δm[i=no_slack_gen], -gp_cc[i] + z_ϵ*t_p[i] <= -generators[i].Pgmin)
    else
        @expression(m_cc, t_p[i=no_slack_gen], 0)
    end


    ## Objective 

    # Quadratic Cost Component
    @variable(m_cc, r_sched >= 0)
    @constraint(m_cc, vcat(r_sched, 0.5, C_rt'*gp_cc) in RotatedSecondOrderCone())
    
    # Linear Cost Component
    @expression(m_cc, linear_cost, sum(generators[i].pi2 * gp_cc[i] + generators[i].pi3 for i in 1:n_generators))

    # Uncertain Component
    if alpha_mod == "split"
        @variable(m_cc, r_bal[i=1:n_generators] >= 0)
        for i in 1:n_generators
            @constraint(m_cc, vec(vcat(r_bal[i], 0.5, (α[i,:]'*Σ_rt)')) in RotatedSecondOrderCone())
        end
        @expression(m_cc, quad_cost, r_sched + sum(r_bal[i]*generators[i].pi1 for i in 1:n_generators))
    else
        @variable(m_cc, r_bal >= 0)
        @constraint(m_cc, vec(vcat(r_bal, 0.5, C_rt'*α*s)) in RotatedSecondOrderCone())
        @expression(m_cc, quad_cost, r_sched + r_bal)
    end


    # Variance Penalty
    @variable(m_cc, zv[1:length(PQ)] >= 0)
    @constraint(m_cc, [i=1:length(PQ)], zv[i] >= t_v[i]^2)
    @variable(m_cc, zq[1:length(PV)] >= 0)
    @constraint(m_cc, [i=1:length(PV)], zq[i] >= t_qPV[i]^2)
    @variable(m_cc, zp[no_slack_gen] >= 0)
    @constraint(m_cc, [i=no_slack_gen], zp[i] >= t_p[i]^2)
    @variable(m_cc, zfp_ij[i=1:n_lines] >= 0)
    @variable(m_cc, zfp_ji[i=1:n_lines] >= 0)
    @constraint(m_cc, [i=no_slack_gen], zfp_ij[i] >= t_f_pij[i]^2)
    @constraint(m_cc, [i=no_slack_gen], zfp_ji[i] >= t_f_pji[i]^2)
    @variable(m_cc, zfq_ij[i=1:n_lines] >= 0)
    @variable(m_cc, zfq_ji[i=1:n_lines] >= 0)
    @constraint(m_cc, [i=no_slack_gen], zfq_ij[i] >= t_f_qij[i]^2)
    @constraint(m_cc, [i=no_slack_gen], zfq_ji[i] >= t_f_qji[i]^2)

    @expression(m_cc, var_penalty,  var_penalty["v"] * sum(zv)
                                    + var_penalty["q_G"] * sum(zq)
                                    + var_penalty["p_G"] * sum(zp)
                                    + var_penalty["f_p"] * (sum(zfp_ij) + sum(zfp_ji))
                                    + var_penalty["f_q"] * (sum(zfq_ij) + sum(zfq_ji)))


    # Objective Funcion`
    @objective(m_cc, Min, quad_cost + linear_cost + var_penalty)

    out_data = Dict(
        # "VSP" => VSP,
        # "VSQ" => VSQ,
        # "θSP_PQ" => θSP_PQ,
        # "θSP_PV" => θSP_PV,
        # "θSQ_PQ" => θSQ_PQ,
        # "θSQ_PV" => θSQ_PV,
        "PQ" => PQ, 
        "PV" => PV,
        "REF" => REF,
        # "A" => A,
        # "C" => C
        )

    return m_cc, out_data

end