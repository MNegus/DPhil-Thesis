function SolStruct = NormalModesTemporalForm(ts, ...
    SolStruct_d_form, alpha, epsilon, DELTA_T)

    % Loads solutions from 
    t_vals_d_form = SolStruct_d_form.t_vals;
    d_vals_d_form = SolStruct_d_form.d_vals;
    as_d_form = SolStruct_d_form.a_vals;
    bs_d_form = SolStruct_d_form.b_vals;
    a_ts_d_form = SolStruct_d_form.a_t_vals;
    d_ts_d_form = SolStruct_d_form.d_t_vals;
    Js_d_form = SolStruct_d_form.J_vals;
    omegas = SolStruct_d_form.omegas;
    E_outers_d_form = SolStruct_d_form.E_outer_vals;

    % Interpolates solutions for time dependent quantities
    SolStruct.ds = interp1(t_vals_d_form, d_vals_d_form, ts);
    SolStruct.as = interp1(t_vals_d_form, as_d_form, ts);
    SolStruct.bs = interp1(t_vals_d_form, bs_d_form, ts);
    SolStruct.a_ts = interp1(t_vals_d_form, a_ts_d_form, ts);
    SolStruct.d_ts = interp1(t_vals_d_form, d_ts_d_form, ts);
    SolStruct.Js = interp1(t_vals_d_form, Js_d_form, ts);
    
    % Determines solution for a_tts by numerical differentiation
    SolStruct.a_tts = zeros(size(SolStruct.a_ts));
%     SolStruct.a_tts(2 : end, :) = diff(SolStruct.a_ts, 1, 1) / DELTA_T;
    SolStruct.a_tts(2 : end - 1, :) = (SolStruct.a_ts(3 : end, :) - SolStruct.a_ts(1 : end - 2, :)) / (2 * DELTA_T);
    SolStruct.a_tts(end, :) = SolStruct.a_tts(end - 1, :);
    
    % Determines q and q_ts using governing equations
    SolStruct.qs = -(alpha / epsilon^2) * (SolStruct.a_ts - omegas' .* SolStruct.bs);
    SolStruct.q_ts = -(alpha / epsilon^2) * (SolStruct.a_tts + omegas' .* SolStruct.as); 

    % Determine energy solutions
    SolStruct.E_outers = interp1(t_vals_d_form, E_outers_d_form, ts);
    E_jets = zeros(size(ts));
    E_jets(2 : end) ...
        = 4 * epsilon^2 * cumtrapz(ts(2 : end), ...
        SolStruct.d_ts(2 : end).^3 .* SolStruct.Js(2 : end));
    SolStruct.E_jets = E_jets;
end

