function [epsilon, k, q, omega] = substrateparameters()
%QUADRATICPARAMETERS Set consistent parameters for moving substrate.
%   The substrate parameters are defined such that we have for the flat
%   substrate,
%   w(x, t) = w(t) = q * (t.^2 + 1 - cos(omega * t)),
%   and for the curved substrate,
%   w(x, t) = w(t) (1 - k^2 * x.^2).
%   
%   Assuming that tmax = 1, the turnover point (in 2D) at tmax is at 
%   epsilon * dmax = 2 * epsilon.
%   If we want the root of the quadratic to remain fixed at some a *
%   epsilon * dmax, then we must have k = 1 / (2 * a * epsilon);

    % Small time parameter
    epsilon = 0.1;
    
    % Oscillation magnitude
    q = 0.1;
    
    % Oscillation frequency
	omega = 3.0;
    
    % Curved substrate root parameter.
%     k = 2 / (3 * epsilon); % a = 3 / 4
    k = 1 / epsilon; % a = 1/2
%     k = 3 / (4 * epsilon); % a = 2 / 3
end

