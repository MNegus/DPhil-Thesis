function Fs = AnalyticalForceAxi(ts, ws, w_ts, w_tts, epsilon)
%%AnalyticalForceAxi
% Determines the analytical force, Fs, in the axisymmetric flat substrate 
% case

    %% Time-dependent quantities
    ds = sqrt(3 * (ts - ws));
    d_ts = sqrt(3) * (1 - w_ts) ./ (2 * sqrt(ts - ws));
    d_tts = -sqrt(3) * ((1 - w_ts).^2 + 2 * (ts - ws) .* w_tts) ./ (4 * (ts - ws).^(3/2));
%     Js = 2 * ds.^3 / (9 * pi);
    Js = 2 * (ts - ws).^(3/2) / (sqrt(3) * pi);
    cs = cSolution(ds, Js, epsilon);
    
    %% Outer force contribution
    FsOuter = (8 * epsilon / 9) * ds.^3 .* (ds .* d_tts + 4 * d_ts.^2);
    
    %% Inner force contribution
    FsInner = (8 * epsilon^4 * d_ts.^2 .* Js.^2 .* cs / pi) ...
        .* (pi * ds ./ (epsilon^2 * Js) - cs.^2 / 3 - 2 * cs - 2 * log(cs) +1);
    
    %% Overlap force contribution
    Fsoverlap = 16 * sqrt(2) * epsilon * ds.^3 .* d_ts.^2 / 9;
    
    %% Composite force
    Fs = FsOuter + FsInner - Fsoverlap;

end

