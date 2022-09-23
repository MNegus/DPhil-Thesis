%% SaveAnalyticalSolutions.m
% Saves all of the analytical solutions used in the plate chapter
clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Forces/");

% General parameters
epsilon = 1;
tMax = 1/12;

% Make directories
mkdir("AnalyticalSolutions");

%% Save alpha varying functions
fineDir = "AnalyticalSolutions/Fine_ALPHA_varying";
mkdir(fineDir);

ALPHAS = 10.^linspace(-2, log10(300), 500);
BETA = 0;
GAMMA = 0;
ALPHATildes = (45 / (8 * epsilon)) * ALPHAS;


for ALPHA = ALPHAS
    ALPHA
    ALPHATilde = (45 / (8 * epsilon)) * ALPHA;

    % Solve asymptotic solution
    xiMax = tMax / ALPHATilde^(2/3);
    DELTA_XI = xiMax / (1e2 - 1);
    xis = 0 : DELTA_XI : xiMax + DELTA_XI;
    fsZeroFun = @(fs) fs + fs.^(2/5) / 3 - xis;
    fsGuess = xis;
    fs = real(fsolve(fsZeroFun, fsGuess));
    
    fPrimesZeroFun = @(fPrimes) fPrimes + 2 * fPrimes ./ (15 * fs(2 : end).^(3/5)) - 1;
    fPrimesGuess = ones(size(xis));
    fPrimes = zeros(size(xis));
    fPrimes(2 : end) = real(fsolve(fPrimesZeroFun, fPrimesGuess(2 : end)));

    % Asymptotic solution
    ts = ALPHATilde^(2/3) * xis;
    ws = ALPHATilde^(2/3) * fs;
    w_ts = fPrimes;
    w_tts = diff(w_ts) ./ diff(ts);

    % Shortens arrays to make same length
    ts = ts(1 : end - 1);
    ws = ws(1 : end - 1);
    w_ts = w_ts(1 : end - 1);
    
    % Find substrate functions
    SubstrateFunctions = platesubstratefunctions(ts, ...
        ws, w_ts, w_tts, epsilon);
    
    % Turnover points
    ds = SubstrateFunctions.d(ts);
    d_ts = SubstrateFunctions.d_t(ts);
    d_ts(1) = 0;
    
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
    fileName = append(fineDir, "/ALPHA_Asy_", num2str(ALPHA), ".mat");
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
    BETATilde = (45 / (8 * epsilon)) * BETA;
    
    % Save analytical solution
    ts = linspace(0, tMax, 1e2);
    ws = (45 * sqrt(3) / (2 * BETATilde)) * ts.^(3/2);
    w_ts = (3 / 2) * (45 * sqrt(3) / (2 * BETATilde)) * ts.^(1/2);
    w_tts = (1 / 2) * (3 / 2) * (45 * sqrt(3) / (2 * BETATilde)) * ts.^(-1/2);
    w_tts(1) = 0;

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
    fileName = append(fineDir, "/BETA_Asy_", num2str(BETA), ".mat");
    save(fileName, "SolStruct");
    
    % Clear workspace variable
    clear("SolStruct");

end

%% Save gamma varying functions
fineDir = "AnalyticalSolutions/Fine_GAMMA_varying";
mkdir(fineDir);

ALPHA = 2;
ALPHATilde = (45 / (8 * epsilon)) * ALPHA;
BETA = 0;
GAMMAS = 10.^linspace(-1, 7, 500);

% for ALPHA_str = ALPHA_strs
for GAMMA = GAMMAS
    GAMMA
    GAMMATilde = (45 / (8 * epsilon)) * GAMMA;

    % Save asymptotic solutions
    xiMax = tMax * sqrt(GAMMATilde / ALPHATilde);
    DELTA_XI = xiMax / (1e2 - 1);
    xis = 0 : DELTA_XI : xiMax + 2 * DELTA_XI;
    gs = zeros(size(xis));
    for xiIdx = 2 : length(xis)
        xi = xis(xiIdx);
        taus = xis(1 : xiIdx);
        gs(xiIdx) = trapz(taus, taus.^0.5 .* sin(xi - taus));
    end
    ts = sqrt(ALPHATilde / GAMMATilde) * xis;
    ws = (135 * sqrt(3) / 4) * (ALPHATilde / GAMMATilde^5)^(1/4) * gs;
    w_ts = diff(ws) ./ diff(ts);
    w_tts = diff(w_ts) ./ diff(ts(1 : end - 1));
    ts = ts(1 : end - 2);
    ws = ws(1 : end - 2);
    w_ts = w_ts(1 : end - 1);

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
    fileName = append(fineDir, "/GAMMA_Asy_", num2str(GAMMA), ".mat");
    save(fileName, "SolStruct");
    
    % Clear workspace variable
    clear("SolStruct");

end
