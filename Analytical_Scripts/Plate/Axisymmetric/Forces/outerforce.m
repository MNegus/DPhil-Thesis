function Fs = outerforce(ts, ws, w_ts, w_tts, epsilon)
%%outerforce
%   Returns the outer solution for the force on the substrate in the
%   axisymmetric-dimensional, plate impact case.

    Fs = 2 * sqrt(3) * epsilon ...
        * (3 * (1 - w_ts).^2 - 2 * w_tts .* (ts - ws)) .* sqrt(ts - ws);
end