function [p, d, d_t, J] = w_dependents(xs, t, w_fun, ...
    w_t_fun, w_tt_fun, w_x_fun, pressure_type, epsilon)
    
    %% Determine d and d_t
    [d, d_t] = TurnoverPointFD(t, w_fun, w_t_fun, w_x_fun, epsilon);

    % Finds d_idx such that x(d_idx) < epsilon * d but x(d_idx) >= epsilon * d
    d_idx = sum(xs < epsilon * d);

    %% Determine m and its derivatives
    if (d_idx < 2)
        m_t_fun = @(s) zeros(size(s));
        m_tt_fun = @(s) zeros(size(s));
    else
        s_vals = xs(1 : d_idx) / epsilon;

        m_t = cumtrapz(s_vals, w_t_fun(epsilon * s_vals), 1);
        m_tt = cumtrapz(s_vals, w_tt_fun(epsilon * s_vals), 1);

        m_t_fun = @(s) interp1(s_vals, m_t, s, 'linear', 'extrap');
        m_tt_fun = @(s) interp1(s_vals, m_tt, s, 'linear', 'extrap');
    end
    
    %% Determine time-dependent quantities
    [A, B, C, J] = TimeDependentsFD(d, d_t, w_t_fun, w_tt_fun, m_t_fun, epsilon);
    
    %% Determine pressure at current timestep
    if pressure_type == "outer"
        p = outer_pressure(xs, m_tt_fun, w_tt_fun, d, A, epsilon);
    elseif pressure_type == "composite"
        p = composite_pressure(xs, t, d, d_t, A, C, J, m_tt_fun, w_tt_fun, epsilon);
    else
        error("Invalid pressure_type");
    end
    
    % Restrict large values of p
    p(p > 1e4) = 0;

end