function TimeDependents = timedependents(ts, SubstrateCoefficients)
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
    ws = SubstrateCoefficients.ws;
    w_ts = SubstrateCoefficients.w_ts;
    w_tts = SubstrateCoefficients.w_tts;
    
    %% Turnover curves
    ds = sqrt(3 * (ts - ws));
    d_ts = sqrt(3) * (1 - w_ts) ./ (2 * sqrt(ts - ws));
    d_tts = - sqrt(3) * ((1 - w_ts).^2 + 2 * (ts - ws) .* w_tts) ./ (4 * (ts - ws).^(3/2));
    
    %% Jet thickness
    Js = 2 * ds.^3 / (9 * pi);
    
    %% Outer-outer pressure coefficient
    As = (4 * ds.^3 / 9) .* (4 * d_ts.^2 + ds .* d_tts);
    
    %% Fill in the structure
    TimeDependents.ds = ds;
    TimeDependents.d_ts = d_ts;
    TimeDependents.d_tts = d_tts;
    TimeDependents.Js = Js;
    TimeDependents.As = As;
    
end