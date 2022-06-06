function [ws, w_ts, w_tts] = ImposedPlate(ts, a)
%ImposedPlate
% Script to return the time-dependence of the imposed plate motion used in
% the validation for the DNS.

    %% Parameters
    omega = 12; % Frequency of the cosine term
    
    %% Solution for ws
    ws = a * (1 - cos(omega *  ts) + ts.^2);
    
    %% Solution for w_ts (first time derivative)
    w_ts = a * (2 * ts + omega * sin(omega * ts));
    
    %% Solution for w_tts (second time derivative)
    w_tts = a * (2 + omega^2 * cos(omega * ts));

end