function [ds, as, a_ts, a_tts, q_ts, d_ts, Js] = NormalModesTemporalForm(ts, ...
    t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, ...
    d_ts_d_form, Js_d_form, omegas, alpha, epsilon, DELTA_T)

    % Interpolates solutions for ds, as and a_ts
    ds = interp1(t_vals_d_form, d_vals_d_form, ts);
    as = interp1(t_vals_d_form, as_d_form, ts);
    a_ts = interp1(t_vals_d_form, a_ts_d_form, ts);
    d_ts = interp1(t_vals_d_form, d_ts_d_form, ts);
    Js = interp1(t_vals_d_form, Js_d_form, ts);
    
    % Determines solution for a_tts by numerical differentiation
    a_tts = zeros(size(a_ts));
    a_tts(2 : end, :) = diff(a_ts, 1, 1) / DELTA_T;
    
    % Determines q_ts using governing equations
    q_ts = -(alpha / epsilon)^2 * (a_tts + omegas' .* as); 
end

