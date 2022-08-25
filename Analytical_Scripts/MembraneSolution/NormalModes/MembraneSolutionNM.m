function [ws, w_ts, ps] = MembraneSolutionNM(xs, as, a_ts, q_ts, d, L, N, epsilon)
%MEMBRANESOLUTIONNM Outputs the solution for the membrane with normal modes
%   At a particular time, for N normal modes, outputs the membrane
%   displacement, ws, its time deriative w_ts, and the pressure ps, given
%   the length N vectors as, a_ts, q_ts, and the turnover point d.
    
    %% Save lambdas
    lambdas = pi * (2 * (1 : N) - 1) / (2 * L);
    
    %% Determine the solutions as a function of xs
    % Membrane displacement
    ws = sum(as .* cos(xs * lambdas), 2) / sqrt(L);

    % Membrane time derivative
    w_ts = sum(a_ts .* cos(xs * lambdas), 2) / sqrt(L);

    % Pressure
    ps = sum(-q_ts .* cos(xs * lambdas), 2) / sqrt(L);
    ps(xs >= epsilon * d) = 0; % Cut off pressure for x > epsilon * d
    

end