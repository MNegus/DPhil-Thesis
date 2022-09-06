function SaveNMSolution(data_dir, alpha, beta, gamma, epsilon, N, L, tmax, delta_t)
    
    %% Derived parameters
    ts = 0 : delta_t : tmax;

    %% d dependent parameters
    d_max = 2 * sqrt(tmax);
    delta_d = delta_t;

    %% Solves the ode in d-form
    NM_Struct_d_form ...
        = NormalModesODE(alpha, beta, gamma, epsilon, delta_d, d_max, N, L);

    %% Converts solution to t-form
    NM_Struct = NormalModesTemporalForm(ts, ...
        NM_Struct_d_form, alpha, delta_t);
    
    %% Saves matrices in data_dir
    SolStruct.N = N;
    SolStruct.ts = NM_Struct.ts;
    SolStruct.ds = NM_Struct.ds;
    SolStruct.as = NM_Struct.as;
    SolStruct.bs = NM_Struct.bs;
    SolStruct.a_ts = NM_Struct.a_ts;
    SolStruct.a_tts = NM_Struct.a_tts;
    NM_Struct.qs = NM_Struct.qs;
    SolStruct.q_ts = NM_Struct.q_ts;
    SolStruct.d_ts = NM_Struct.d_ts;
    SolStruct.Js = NM_Struct.Js;
    SolStruct.E_outers = NM_Struct.E_outers;
    SolStruct.E_jets = NM_Struct.E_jets;

    save(sprintf("%s/SolStruct.mat", data_dir), 'SolStruct');
end