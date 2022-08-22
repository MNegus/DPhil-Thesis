%% SaveAnalyticalSolutions.m
% Saves all of the analytical solutions used in the plate chapter

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Forces/");

% General parameters
epsilon = 1;
tMax = 1/3;

% Make directories
mkdir("AnalyticalSolutions");

%% Save alpha varying functions
fineDir = "AnalyticalSolutions/Fine_ALPHA_varying";
mkdir(fineDir);

% ALPHAS = linspace(1, 200, 100);
ALPHAS = 10.^linspace(-2, log10(300), 500);
BETA = 0;
GAMMA = 0;

% for ALPHA_str = ALPHA_strs
for ALPHA = ALPHAS
    ALPHA
    % Save plate solution 
    [ts, ws, w_ts, w_tts] ...
        = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

    % Find substrate functions
    SubstrateFunctions = platesubstratefunctions(ts, ...
        ws, w_ts, w_tts, epsilon);
    
    % Turnover points
    ds = SubstrateFunctions.d(ts);
    d_ts = SubstrateFunctions.d_t(ts);
    
    % Force solution
    [FsComp, FsOuter, ~] ...
        = substrateforce(ts, SubstrateFunctions);
    
    % Jet thickness
    Js = SubstrateFunctions.J(ts);

    % Energies
    [EOuters, EJets] = dropletenergy(ts, SubstrateFunctions);

    % Create solutions struct
    SolStruct.ts = ts;
    SolStruct.ws = ws;
    SolStruct.w_ts = w_ts;
    SolStruct.w_tts = w_tts;
    SolStruct.ds = ds;
    SolStruct.d_ts = d_ts;
    SolStruct.FsComp = FsComp;
    SolStruct.Js = Js;
    SolStruct.EOuters = EOuters;
    SolStruct.EJets = EJets;
    SolStruct.SubstrateFunctions = SubstrateFunctions;

    % Save solutions struct
    fileName = append(fineDir, "/ALPHA_", num2str(ALPHA), ".mat");
    save(fileName, "SolStruct");
    
    % Clear workspace variable
    clear("SolStruct");

end


%% Save beta varying functions
fineDir = "AnalyticalSolutions/Fine_BETA_varying";
mkdir(fineDir);

ALPHA = 2;
GAMMA = 100;
BETA_crit = 2 * sqrt(ALPHA * GAMMA);
BETAS = BETA_crit * 10.^linspace(-3, 3, 500);

% for ALPHA_str = ALPHA_strs
for BETA = BETAS
    BETA
    % Save plate solution 
    [ts, ws, w_ts, w_tts] ...
        = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

    % Find substrate functions
    SubstrateFunctions = platesubstratefunctions(ts, ...
        ws, w_ts, w_tts, epsilon);
    
    % Turnover points
    ds = SubstrateFunctions.d(ts);
    d_ts = SubstrateFunctions.d_t(ts);
    
    % Force solution
    [FsComp, FsOuter, ~] ...
        = substrateforce(ts, SubstrateFunctions);
    
    % Jet thickness
    Js = SubstrateFunctions.J(ts);

    % Energies
    [EOuters, EJets] = dropletenergy(ts, SubstrateFunctions);

    % Create solutions struct
    SolStruct.ts = ts;
    SolStruct.ws = ws;
    SolStruct.w_ts = w_ts;
    SolStruct.w_tts = w_tts;
    SolStruct.ds = ds;
    SolStruct.d_ts = d_ts;
    SolStruct.FsComp = FsComp;
    SolStruct.Js = Js;
    SolStruct.EOuters = EOuters;
    SolStruct.EJets = EJets;
    SolStruct.SubstrateFunctions = SubstrateFunctions;

    % Save solutions struct
    fileName = append(fineDir, "/BETA_", num2str(BETA), ".mat");
    save(fileName, "SolStruct");
    
    % Clear workspace variable
    clear("SolStruct");

end

%% Save gamma varying functions
fineDir = "AnalyticalSolutions/Fine_GAMMA_varying";
mkdir(fineDir);

ALPHA = 2;
BETA = 0;
GAMMAS = 10.^linspace(-1, 7, 500);

% for ALPHA_str = ALPHA_strs
for GAMMA = GAMMAS
    GAMMA
    % Save plate solution 
    [ts, ws, w_ts, w_tts] ...
        = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

    % Find substrate functions
    SubstrateFunctions = platesubstratefunctions(ts, ...
        ws, w_ts, w_tts, epsilon);
    
    % Turnover points
    ds = SubstrateFunctions.d(ts);
    d_ts = SubstrateFunctions.d_t(ts);
    
    % Force solution
    [FsComp, FsOuter, ~] ...
        = substrateforce(ts, SubstrateFunctions);
    
    % Jet thickness
    Js = SubstrateFunctions.J(ts);

    % Energies
    [EOuters, EJets] = dropletenergy(ts, SubstrateFunctions);

    % Create solutions struct
    SolStruct.ts = ts;
    SolStruct.ws = ws;
    SolStruct.w_ts = w_ts;
    SolStruct.w_tts = w_tts;
    SolStruct.ds = ds;
    SolStruct.d_ts = d_ts;
    SolStruct.FsComp = FsComp;
    SolStruct.Js = Js;
    SolStruct.EOuters = EOuters;
    SolStruct.EJets = EJets;
    SolStruct.SubstrateFunctions = SubstrateFunctions;

    % Save solutions struct
    fileName = append(fineDir, "/GAMMA_", num2str(GAMMA), ".mat");
    save(fileName, "SolStruct");
    
    % Clear workspace variable
    clear("SolStruct");

end
