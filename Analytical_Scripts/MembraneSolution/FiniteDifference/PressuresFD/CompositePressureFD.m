function ps = CompositePressureFD(xs, t, d, d_t, A, C, J, m_tt_fun, w_tt_fun, epsilon)

    if (t == 0)
        ps = zeros(size(xs));
    else
        % Calculate outer pressure
        outer_ps = OuterPressureFD(xs, m_tt_fun, w_tt_fun, d, A, epsilon);
        
        % Calculate overlap pressure
        overlap_ps = OverlapPressureFD(xs, d, C, epsilon);

        % Calculate inner pressure
        inner_ps = InnerPressureFD(xs, d, d_t, J, epsilon);

        % Calculate composite pressure
        ps = outer_ps + inner_ps - overlap_ps;
    end
    
end

