function SaveValidatedNMSolution(data_dir, alpha, beta, gamma, epsilon, L, tmax, delta_t)

    [N, delta_d, ts, ds, as, a_ts, a_tts, q_ts, d_ts, Js] ...
        = ValidatedNMSolution(alpha, beta, gamma, epsilon, L, tmax, delta_t);
    
    %% Saves matrices in data_dir
    SolStruct.N = N;
    SolStruct.delta_d = delta_d;
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