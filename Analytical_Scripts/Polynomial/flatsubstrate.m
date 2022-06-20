function [as, a_ts, a_tts, bs, b_ts, b_tts] ...
    = flatsubstrate(ts, q, omega)
%FLATSUBSTRATE Returns substrate coefficient functions for a flat
%substrate.
%   We assume a quadratic substrate of the form
%   w(x, t) = w(t),
%   where
%   w(t) = q * (t + 1 - cos(omega * t)),
%   we return the substrate coefficients such that
%   w(x, t) = a(t).

    % a coefficients
    as = q * (ts.^2 + 1 - cos(omega * ts));
    a_ts = q * (2 * ts + omega * sin(omega * ts));
    a_tts = q * (2 + omega^2 * cos(omega * ts));

    % b coefficients
    bs = zeros(size(as));
    b_ts = zeros(size(as));
    b_tts = zeros(size(as));
end
