function TimeDependents = timedependents(ts, ws, w_ts, w_tts)
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

    % Define structure array
    

    %% Turnover points
    ds = 2 * sqrt(ts - ws);
    d_ts = (1 - w_ts) ./ sqrt(ts - ws);
    
    %% As, Bs and Cs coefficients
    As = ds .* d_ts .* (1 - w_ts) - 0.5 * ds.^2 .* w_tts;
    Bs = w_ts;
    Cs = ds .* d_ts .* (1 - Bs);
    
    %% Jet thickness
    Js = pi * (1 - Bs).^2 .* ds ./ (8 * d_ts.^2);
    
    %% Fill in the structure
    TimeDependents.ds = ds;
    TimeDependents.d_ts = d_ts;
    TimeDependents.As = As;
    TimeDependents.Bs = Bs;
    TimeDependents.Cs = Cs;
    TimeDependents.Js = Js;
    
end