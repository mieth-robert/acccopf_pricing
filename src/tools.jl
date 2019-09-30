raw_index(v::MOI.VariableIndex) = v.value
raw_index(c::NonlinearConstraintIndex) = c.value


function extract_jacobian(m, nonlinear_constraints, variables)
    all_vars = all_variables(m)
    n_vars = length(all_vars)
    const_types = list_of_constraint_types(m)
    n_allconstr = sum(num_constraints(m, ct[1], ct[2]) for ct in const_types)

    d = NLPEvaluator(m)
    requested_features = [:Jac, :Grad] 
    MOI.initialize(d, requested_features)
    jac_struct = MOI.jacobian_structure(d)

    current_results = value.(all_vars)
    V = zeros(length(jac_struct))

    MOI.eval_constraint_jacobian(d, V, current_results)

    I = [js[1] for js in jac_struct] # rows
    J = [js[2] for js in jac_struct] # cols

    # jac = sparse(I, J, V, n_allconstr, n_vars)
    jac = sparse(I, J, V)

    grad = zeros(n_allconstr)
    MOI.eval_objective_gradient(d, grad, current_results)

    rows = raw_index.(index.(nonlinear_constraints))
    cols = raw_index.(index.(variables))

    return  grad[rows], jac[rows, cols]
end
