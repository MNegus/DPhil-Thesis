function Fs = overlapforce(ts, ws, w_ts, epsilon)
%%innerforce
%   Returns the overlap solution for the force on the substrate in the
%   two-dimensional, plate impact case.
    
    Fs = 4 * sqrt(6) * epsilon * (1 - w_ts).^2 .* sqrt(ts - ws);
end