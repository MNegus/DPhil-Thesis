function SaveFDSolution(parent_dir, ...
    ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, ...
    pressure_type)

    %% Saves data in sub-directory depending on pressure type
    fd_data_dir = sprintf("%s/%s", parent_dir, pressure_type);

    %% Derived parameters
    DELTA_X = L / (N_MEMBRANE - 1); 
    M = N_MEMBRANE - 1;
    xs = (0 : DELTA_X : L - DELTA_X)';
    T_VALS = 0 : DELTA_T : T_MAX;
    
    if (GAMMA == 0)
        Cpressure = 4 * DELTA_X * DELTA_X / BETA;
    else
        Cpressure = 4 * DELTA_X^4 / GAMMA;
    end
    

    %% Matrix definitions
    [A_mat, B_mat] = MatrixDefinitions(ALPHA, BETA, GAMMA, M, DELTA_X, DELTA_T);

    %% Initialise arrays
    w_previous = zeros(size(xs));

    % Save initial array
    w_next = w_previous;
    save(sprintf("%s/w_%d.mat", fd_data_dir, 0), 'w_next');

    % Save array at k = 1
    w = InitialiseMembrane(w_previous, A_mat, B_mat);
    w_next = w;
    save(sprintf("%s/w_%d.mat", fd_data_dir, 1), 'w_next');

    % Initialise remaining arrays
    w_t = zeros(size(xs));
    save(sprintf("%s/w_t_%d.mat", fd_data_dir, 0), 'w_t');
    p = zeros(size(xs));
    save(sprintf("%s/p_%d.mat", fd_data_dir, 0), 'p');
    p_previous = zeros(size(xs));
    ds = zeros(size(T_VALS));
    d_ts = zeros(size(T_VALS));
    Js = zeros(size(T_VALS));
    pMaxs = zeros(size(T_VALS));
    p0s = zeros(size(T_VALS));
    
    %% Loops over time
    for k = 1 : length(T_VALS) - 1
        %% Updates time
        t = T_VALS(k);
        t

        %% Composite timestep
        [w_next, p, w_t, d, d_t, J, pMax] = MembraneTimestep(...
            xs, t, w, w_previous, p_previous, pressure_type, ...
            EPSILON, DELTA_T, DELTA_X, M, Cpressure, A_mat, B_mat);
        
%         plot(xs, w);
%         drawnow;
        
        ds(k) = d;
        d_ts(k) = d_t;
        Js(k) = J;
        pMaxs(k) = pMax;
        p0s(k) = p(1);

        % Saves arrays
        save(sprintf("%s/w_%d.mat", fd_data_dir, k + 1), 'w_next');
        save(sprintf("%s/w_t_%d.mat", fd_data_dir, k), 'w_t');
        save(sprintf("%s/p_%d.mat", fd_data_dir, k), 'p');
        
        % Swaps ws
        temp = w_previous;
        w_previous = w;
        w = w_next;
        w_next = temp;
        
        % Swaps pressure
        p_previous = p;

    end
    
    %% Saves ds solutions
    save(sprintf("%s/ds.mat", fd_data_dir), 'ds');
    save(sprintf("%s/d_ts.mat", fd_data_dir), 'd_ts');
    save(sprintf("%s/Js.mat", fd_data_dir), 'Js');
    save(sprintf("%s/pMaxs.mat", fd_data_dir), 'pMaxs');
    save(sprintf("%s/p0s.mat", fd_data_dir), 'p0s');
    
end