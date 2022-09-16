 function ps = OuterPressureFD(xs, m_tt_fun, w_tt_fun, d, A, epsilon)
    ps = zeros(size(xs));
    
    if (d == 0)
        return;
    end
    
    xhats = xs / epsilon;

    %% Finds idx such that xhat = d
    idx = sum(xhats < d);
    
    %% Determines s values
    s_vals = xhats(1 : idx)';
    
    %% Full matrix method
    xhat_dependents_alt = (s_vals .* m_tt_fun(s_vals))';
    
    full_integrands = trapz_matrix(s_vals, xhats(1 : idx), d, m_tt_fun, w_tt_fun, epsilon);
    integrals = trapz(s_vals, full_integrands, 2);
    xhat_dependents_alt = xhat_dependents_alt - integrals;
    
    %% Returns p values
    ps(1 : idx) = (A + xhat_dependents_alt) ./ sqrt(epsilon^2 * d^2 - xs(1 : idx).^2);
    
    %% Function definitions
    function mat = trapz_matrix(s_vals, xhats, d, m_tt_fun, w_tt_fun, epsilon)
        
        % Sets the bulk nodes
        mat = (2 / pi) * (sqrt(d^2 - s_vals.^2) .* (s_vals .* m_tt_fun(s_vals) ...
            - xhats .* m_tt_fun(xhats))) ./ (s_vals.^2 - (xhats).^2);
        
        % Set the diagonal to zero (as it would be NaN)
        mat(1 : 1 + size(mat,1) : end) = 0;
        
        % Sets the diagonals
        diagonal = (2 / pi) * sqrt(d^2 - xhats.^2) .* (m_tt_fun(xhats) + xhats .* w_tt_fun(epsilon * xhats)) ./ (2 * xhats);
        mat = mat + diag(diagonal);
        mat(1, 1) = (2 / pi) * d * w_tt_fun(0);
    end

end