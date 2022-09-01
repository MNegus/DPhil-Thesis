function SaveNMSolution(data_dir, alpha, beta, gamma, epsilon, N, L, tmax, delta_t)
    
    %% Derived parameters
    ts = 0 : delta_t : tmax;

    %% d dependent parameters
    d_max = 2 * sqrt(tmax);
    delta_d = delta_t;

    %% Solves the ode in d-form
    [t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, kvals] ...
        = NormalModesODE(alpha, beta, gamma, epsilon, delta_d, d_max, N, L);

    %% Converts solution to t-form
    [ds, as, a_ts, a_tts, q_ts, d_ts, Js] = NormalModesTemporalForm(ts, ...
        t_vals_d_form, d_vals_d_form, as_d_form, a_ts_d_form, kvals, alpha, delta_t);
    
    %% Saves matrices in data_dir
    SolStruct.N = N;
    SolStruct.ts = ts;
    SolStruct.ds = ds;
    SolStruct.as = as;
    SolStruct.a_ts = a_ts;
    SolStruct.a_tts = a_tts;
    SolStruct.q_ts = q_ts;
    SolStruct.d_ts = d_ts;
    SolStruct.Js = Js;

    save(sprintf("%s/SolStruct.mat", data_dir), 'SolStruct');
end