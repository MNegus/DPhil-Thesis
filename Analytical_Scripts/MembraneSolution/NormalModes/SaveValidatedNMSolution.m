function SaveValidatedNMSolution(data_dir, alpha, beta, gamma, epsilon, L, tmax, delta_t)

    [N, delta_d, ts, NM_Struct] ...
        = ValidatedNMSolution(alpha, beta, gamma, epsilon, L, tmax, delta_t);
    
    %% Saves matrices in data_dir
    SolStruct.N = N;
    SolStruct.delta_d = delta_d;
    SolStruct.ts = ts;
    SolStruct.ds = NM_Struct.ds;
    SolStruct.as = NM_Struct.as;
    SolStruct.bs = NM_Struct.bs;
    SolStruct.a_ts = NM_Struct.a_ts;
    SolStruct.a_tts = NM_Struct.a_tts;
    SolStruct.qs = NM_Struct.qs;
    SolStruct.q_ts = NM_Struct.q_ts;
    SolStruct.d_ts = NM_Struct.d_ts;
    SolStruct.Js = NM_Struct.Js;
    SolStruct.E_outers = NM_Struct.E_outers;
    SolStruct.E_jets = NM_Struct.E_jets;

    save(sprintf("%s/SolStruct.mat", data_dir), 'SolStruct');

end