%%

% Data directory
dns_master_dir = "/home/michael/scratch/DPhil_DNS_Data/";
dns_dirs = ["Stationary_Plate/axi/", "Moving_Plate/ALPHA-2_BETA-0_GAMMA-20/"];

for dir = dns_dirs
    data_dir = sprintf("%s/%s", dns_master_dir, dir);

    % Timesteps 
    timesteps = [13, 20, 50, 80];
    
    % Number of elements per dimension
    N = 2731;
    
    % Original axes limits
    box_width = 2;
    
    % Coarsened axes limits
    xMax = 2;
    yMax = 2;
    
    % Dimension of coarsened
    NCoarse = 256;
    
    % Fine mesh grid
    xsFull = linspace(0, xMax, N);
    ysFull = linspace(0, yMax, N);
    [XFull, YFull] = meshgrid(xsFull, ysFull);
    
    % Coarse mesh grid
    xs = linspace(0, xMax, NCoarse);
    ys = linspace(0, yMax, NCoarse);
    [X, Y] = meshgrid(xs, ys);
    
    %% Loop over timesteps and coarsen data
    for timestep = timesteps
        timestep
        
        % Read fine data
        data = readmatrix(sprintf("%s/raw_data/field_output_%d.txt", data_dir, timestep));
    
        % Save values (transposing as y is x in Basilisk)
        psFull = data(:, 3)';
        psMesh = reshape(psFull, [N, N]);
    
        fsFull = data(:, 4)';
        fsMesh = reshape(fsFull, [N, N]);
        
        usFull = data(:, 5)';
        usMesh = reshape(usFull, [N, N]);
    
        vsFull = data(:, 6)';
        vsMesh = reshape(vsFull, [N, N]);
    
        % Interpolate onto coarse grid
        ps = interp2(XFull, YFull, psMesh, X, Y);
        fs = interp2(XFull, YFull, fsMesh, X, Y);
        us = interp2(XFull, YFull, usMesh, X, Y);
        vs = interp2(XFull, YFull, vsMesh, X, Y);
    
        % Create structure array
        fieldStruct.X = X;
        fieldStruct.Y = Y;
        fieldStruct.ps = ps;
        fieldStruct.fs = fs;
        fieldStruct.us = us;
        fieldStruct.vs = vs;
    
        % Save structure array
        filename = sprintf("%s/coarsenedOutputs/fieldStruct_%d.mat", data_dir, timestep);
        save(filename, "fieldStruct");
    
    end

end