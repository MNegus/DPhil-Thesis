function [as, a_ts, a_tts, bs, b_ts, b_tts] ...
    = quadraticsubstrate(ts, L, q, omega)
%QUADRATICSUBSTRATE Returns substrate coefficient functions for quadratic
%substrate
%   We assume a quadratic substrate of the form
%   w(x, t) = w(t) (1 - x^2 / L^2),
%   where
%   w(t) = q * (t + 1 - cos(omega * t)),
%   we return the substrate coefficients such that
%   w(x, t) = a(t) + b(t) * x^2.

    % a coefficients
    as = q * (ts.^2 + 1 - cos(omega * ts));
    a_ts = q * (2 * ts + omega * sin(omega * ts));
    a_tts = q * (2 + omega^2 * cos(omega * ts));

    % b coefficients
    bs = - as / L^2;
    b_ts = - a_ts / L^2;
    b_tts = - a_tts / L^2;
end

