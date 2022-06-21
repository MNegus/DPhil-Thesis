function [ws, w_ts, w_tts] = flatsubstrate(ts, q, omega)
%FLATSUBSTRATE Returns substrate coefficient functions for a flat
%substrate.
%   We assume a quadratic substrate of the form
%   w(x, t) = w(t),
%   where
%   w(t) = q * (t + 1 - cos(omega * t)).

    % Coefficients
    ws = q * (ts.^2 + 1 - cos(omega * ts));
    w_ts = q * (2 * ts + omega * sin(omega * ts));
    w_tts = q * (2 + omega^2 * cos(omega * ts));
end
