function Fs = AnalyticalForce2D(ts, ws, w_ts, w_tts, epsilon)
%%AnalyticalForce2D
% Determines the analytical force, Fs, in the 2D flat substrate case

    %% Time-dependent quantities
    ds = 2 * sqrt(ts - ws);
    d_ts = (1 - w_ts) ./ sqrt(ts - ws);
    Js = pi * (1 - w_ts).^2 .* ds ./ (8 * d_ts.^2);
    B = w_ts;
    cs = cSolution(ds, Js, epsilon);
    
    %% Outer force contribution
    FsOuter = pi * (ds .* d_ts .* (1 - w_ts) - 0.5 * ds.^2 .* w_tts);
    
    %% Inner force contribution
    FsInner = 8 * epsilon * d_ts.^2 .* Js .* cs / pi;
    
    %% Overlap force contribution
    Fsoverlap = 2 * sqrt(2) * ds .* d_ts .* (1 - B);
    
    %% Composite force
    Fs = FsOuter + FsInner - Fsoverlap;

end