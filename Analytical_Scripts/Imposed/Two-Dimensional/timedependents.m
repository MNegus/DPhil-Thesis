function TimeDependents = timedependents(ts, SubstrateFunctions)
%%timedependents
% Quadratic substrate case:
% Given the times and the substrate coefficients, returns the 
% time-dependent quantities in a structure array:
%   ds: Turnover point
%   d_ts: Time derivative of ds
%   Js: Jet thicknesses
%   As: Coefficient for outer pressure
%   Bs: Coefficient for jet thickness
%   Cs: Coefficient for overlap pressure

    %% Load in substrate coefficients
    aHats = SubstrateCoefficients.aHat(ts);
    aHat_ts = SubstrateCoefficients.aHat_t(ts);
    aHat_tts = SubstrateCoefficients.aHat_tt(ts);
    
    bHats = SubstrateCoefficients.bHat(ts);
    bHat_ts = SubstrateCoefficients.bHat_t(ts);
    bHat_tts = SubstrateCoefficients.bHat_tt(ts);
    
    %% Turnover points
    ds = 2 * sqrt((ts - aHats) ./ (1 + 2 * bHats));
    d_ts = ((1 + 2 * bHats) .* (1 - aHat_ts) - 2 * (ts - aHats) .* bHat_ts) ...
        ./ (sqrt(ts - aHats) .* (1 + 2 * bHats).^(3/2));
    d_tts = -((1 + 2 * bHats) .* d_ts.^2 ...
        + 4 * bHat_ts .* ds .* d_ts ...
        + 2 * aHat_tts ...
        + bHat_tts .* ds.^2) ./ ((1 + 2 * bHats) .* ds);
    
    %% As, Bs and Cs coefficients
    % ENSURE WE COERRECTLY HANDLE AT t = 0 !!!!
    
    % Bs, jet-thickness factor
    Bs = aHat_ts + 0.5 * bHat_ts .* ds.^2;
    
    % Cs, overlap pressure factor
    Cs = ds .* d_ts .* (1 - Bs);
    zeroIdx = find(ts == 0);
    if ~isempty(zeroIdx)
       Cs(zeroIdx) = 2 * (1 - aHat_ts(zeroIdx)); 
    end
    
    % As, outer pressure factor
    As = Cs - 0.5 * aHat_tts .* ds.^2 - (1 / 8) * bHat_tts .* ds.^4;
    
    %% Jet thickness
    Js = pi * (1 - Bs).^2 .* ds ./ (8 * d_ts.^2);
    
    %% Fill in the structure
    TimeDependents.ds = ds;
    TimeDependents.d_ts = d_ts;
    TimeDependents.d_tts = d_tts;
    TimeDependents.As = As;
    TimeDependents.Bs = Bs;
    TimeDependents.Cs = Cs;
    TimeDependents.Js = Js;
    
end