function pMax = MaxPressureFD(t, d, d_t, A, C, J, m_tt_fun, w_tt_fun, epsilon)
    

    if (t == 0)
        pMax = 0;
    else
        % Determine x position of maximum pressure
        xMax = epsilon * d - epsilon^3 * J / pi;

        % Calculate outer pressure
        outer_p = OuterPressureFD(xMax, m_tt_fun, w_tt_fun, d, A, epsilon);
        
        % Calculate overlap pressure
        overlap_p = OverlapPressureFD(xMax, d, C, epsilon);

        % Calculate inner pressure
        inner_p = d_t^2 / (2 * epsilon^2);

        % Calculate composite pressure
        pMax = outer_p + inner_p - overlap_p;
    end
end