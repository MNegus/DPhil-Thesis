%% SaveAnalyticalSolutions.m
% Saves all of the analytical solutions used in the plate chapter

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Forces/");

% General parameters
epsilon = 1;

% Make directories
mkdir("AnalyticalSolutions");

%% Save stationary solutions
tMax = 1 / 3; % Max such that 1 = d(t)
ts = linspace(0, tMax, 1e3);

% Find substrate functions
SubstrateFunctions = platesubstratefunctions(ts, ...
    zeros(size(ts)), zeros(size(ts)), zeros(size(ts)), epsilon);

% Find ds
ds = SubstrateFunctions.d(ts);

% Force solution
[FsComp, FsOuter, ~] ...
    = substrateforce(ts, SubstrateFunctions);

% Jet thickness
Js = SubstrateFunctions.J(ts);

% Create solutions struct
SolStruct.ts = ts;
SolStruct.ws = zeros(size(ts));
SolStruct.w_ts = zeros(size(ts));
SolStruct.w_tts = zeros(size(ts));
SolStruct.ds = ds;
SolStruct.FsComp = FsComp;
SolStruct.Js = Js;
SolStruct.SubstrateFunctions = SubstrateFunctions;

% Save struct
save("AnalyticalSolutions/StationarySol", "SolStruct");

% Clear workspace variable
clear("SolStruct");

%% Save alpha varying functions
mkdir("AnalyticalSolutions/ALPHA_varying");

ALPHA_strs = ["2.0", "5.0", "10.0", "20.0", "100.0"];
BETA = 0;
GAMMA = 0;

for ALPHA_str = ALPHA_strs
    % Load numerical value for ALPHA
    ALPHA = str2double(ALPHA_str);

    % Save plate solution 
    [ts, ws, w_ts, w_tts] ...
        = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

    % Find substrate functions
    SubstrateFunctions = platesubstratefunctions(ts, ...
        ws, w_ts, w_tts, epsilon);
    
    % Find where turnover point reaches 1
    ds = SubstrateFunctions.d(ts);
    tIdxMaxComp = sum(ds <= 1);
    
    % Restrict solutions temporally
    ts = ts(1 : tIdxMaxComp);
    ws = ws(1 : tIdxMaxComp);
    w_ts = w_ts(1 : tIdxMaxComp);
    w_tts = w_tts(1 : tIdxMaxComp);
    ds = ds(1 : tIdxMaxComp);
    
    % Force solution
    [FsComp, FsOuter, ~] ...
        = substrateforce(ts, SubstrateFunctions);
    
    % Jet thickness
    JsAnalytical = SubstrateFunctions.J(ts);

    % Create solutions struct
    SolStruct.ts = ts;
    SolStruct.ws = ws;
    SolStruct.w_ts = w_ts;
    SolStruct.w_tts = w_tts;
    SolStruct.ds = ds;
    SolStruct.FsComp = FsComp;
    SolStruct.Js = Js;
    SolStruct.SubstrateFunctions = SubstrateFunctions;

    % Save solutions struct
    fileName = append("AnalyticalSolutions/ALPHA_varying/ALPHA_", ALPHA_str, ".mat");
    save(fileName, "SolStruct");
    
    % Clear workspace variable
    clear("SolStruct");

end

%% Save beta varying functions
mkdir("AnalyticalSolutions/BETA_varying");

ALPHA = 2;
BETA_strs = ["0.0", "7.07", "28.28", "141.42"];
GAMMA = 100;

for BETA_str = BETA_strs
    % Load numerical value for BETA
    BETA = str2double(BETA_str);

    % Save plate solution 
    [ts, ws, w_ts, w_tts] ...
        = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

    % Find substrate functions
    SubstrateFunctions = platesubstratefunctions(ts, ...
        ws, w_ts, w_tts, epsilon);
    
    % Find where turnover point reaches 1
    ds = SubstrateFunctions.d(ts);
    tIdxMaxComp = sum(ds <= 1);
    
    % Restrict solutions temporally
    ts = ts(1 : tIdxMaxComp);
    ws = ws(1 : tIdxMaxComp);
    w_ts = w_ts(1 : tIdxMaxComp);
    w_tts = w_tts(1 : tIdxMaxComp);
    ds = ds(1 : tIdxMaxComp);
    
    % Force solution
    [FsComp, FsOuter, ~] ...
        = substrateforce(ts, SubstrateFunctions);
    
    % Jet thickness
    JsAnalytical = SubstrateFunctions.J(ts);

    % Create solutions struct
    SolStruct.ts = ts;
    SolStruct.ws = ws;
    SolStruct.w_ts = w_ts;
    SolStruct.w_tts = w_tts;
    SolStruct.ds = ds;
    SolStruct.FsComp = FsComp;
    SolStruct.Js = Js;
    SolStruct.SubstrateFunctions = SubstrateFunctions;

    % Save solutions struct
    fileName = append("AnalyticalSolutions/BETA_varying/BETA_", BETA_str, ".mat");
    save(fileName, "SolStruct");
    
    % Clear workspace variable
    clear("SolStruct");

end

%% Save gamma varying functions
mkdir("AnalyticalSolutions/GAMMA_varying");

ALPHA = 2;
BETA = 0;
GAMMA_strs = ["10.0", "20.0", "100.0", "500.0", "1000.0"];

for GAMMA_str = GAMMA_strs
    % Load numerical value for BETA
    GAMMA = str2double(GAMMA_str);

    % Save plate solution 
    [ts, ws, w_ts, w_tts] ...
        = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

    % Find substrate functions
    SubstrateFunctions = platesubstratefunctions(ts, ...
        ws, w_ts, w_tts, epsilon);
    
    % Find where turnover point reaches 1
    ds = SubstrateFunctions.d(ts);
    tIdxMaxComp = sum(ds <= 1);
    
    % Restrict solutions temporally
    ts = ts(1 : tIdxMaxComp);
    ws = ws(1 : tIdxMaxComp);
    w_ts = w_ts(1 : tIdxMaxComp);
    w_tts = w_tts(1 : tIdxMaxComp);
    ds = ds(1 : tIdxMaxComp);
    
    % Force solution
    [FsComp, FsOuter, ~] ...
        = substrateforce(ts, SubstrateFunctions);
    
    % Jet thickness
    JsAnalytical = SubstrateFunctions.J(ts);

    % Create solutions struct
    SolStruct.ts = ts;
    SolStruct.ws = ws;
    SolStruct.w_ts = w_ts;
    SolStruct.w_tts = w_tts;
    SolStruct.ds = ds;
    SolStruct.FsComp = FsComp;
    SolStruct.Js = Js;
    SolStruct.SubstrateFunctions = SubstrateFunctions;

    % Save solutions struct
    fileName = append("AnalyticalSolutions/GAMMA_varying/GAMMA_", GAMMA_str, ".mat");
    save(fileName, "SolStruct");
    
    % Clear workspace variable
    clear("SolStruct");

end

