function TimeDependents = timedependents(ts, PlatePositions)
%%timedependents
% Given the times and the plate position (ws) and its time derivatives
% (w_ts and w_tts), returns the time-dependent quantities in a structure 
% array:
%   ds: Turnover point
%   d_ts: Time derivative of ds
%   Js: Jet thicknesses
%   As: Coefficient for outer pressure
%   Bs: Coefficient for jet thickness
%   Cs: Coefficient for overlap pressure

    %% Load in plate positions
    ws = PlatePositions.ws;
    w_ts = PlatePositions.w_ts;
    w_tts = PlatePositions.w_tts;
    
    %% Turnover points
    ds = sqrt(3 * (ts - ws));
    d_ts = sqrt(3) * (1 - w_ts) ./ (2 * sqrt(ts - ws));
    d_tts = - sqrt(3) * (2 * (ts - ws).^2 .* w_tts + (1 - w_ts).^2) ...
        ./ (4 * (ts - ws).^(3/2));
    
    %% Jet thickness
    Js = 2 * ds.^3 / (9 * pi);
    
    %% Fill in the structure
    TimeDependents.ds = ds;
    TimeDependents.d_ts = d_ts;
    TimeDependents.d_tts = d_tts;
    TimeDependents.Js = Js;
    
end