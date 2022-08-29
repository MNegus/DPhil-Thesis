function [t_vals, d_vals, a_vals, a_t_vals, omegas] = NormalModesODE(alpha, beta, gamma, epsilon, dd, dmax, N, L)

% Discretised dvals
d_vals = 0 : dd : dmax;

% Solves ODE
ks = pi * (2 * (1 : N)' - 1) / (2 * L);
omegas = epsilon^4 * (beta * ks.^2 + gamma * ks.^4) / alpha;
odefun = @(d, y) full_odefun(d, y, ks, omegas, alpha, epsilon, L, N);
y0 = zeros(2 * N, 1);
opts = odeset('RelTol',1e-4,'AbsTol',1e-4, 'Stats', 'on', 'Maxstep', dd);
[d_vals, y] = ode45(odefun, d_vals, y0, opts);
a_vals = y(:, 1 : N);
bvals = y(:, N + 1 : 2 * N);

% Determines time values
t_vals = zeros(size(d_vals));
for k = 1 : length(t_vals)
    t_vals(k) = d_vals(k)^2 / 4 ...
        + dot(a_vals(k, :), besselj(0, epsilon * d_vals(k) * ks)) / sqrt(L);
end

% Determine a_t vals
a_t_vals = zeros(length(t_vals), N);
for k = 1 : length(t_vals)
    d = d_vals(k);
    gvals = pi * d * besselj(1, epsilon * d * ks) ./ (epsilon * ks);
    Mval = MassMatrix(d, alpha, epsilon, L, N);
    a_t_vals(k, :) = Mval \ (omegas .* bvals(k, :)' + epsilon^2 * gvals / (alpha * sqrt(L)));
end

%% ODE function definition
function dydd = full_odefun(d, y, ks, omegas, alpha, epsilon, L, N)
    d
    
    dydd = zeros(2 * N, 1);
    
    %% Load in values
    as = y(1 : N, 1);
    bs = y(N + 1 : 2 * N, 1);
    
    %% Calculate useful quantities
   
    gs = pi * d * besselj(1, epsilon * d * ks) ./ (epsilon * ks);
    Gammas = besselj(0, epsilon * d * ks);
    Gammas_d = - epsilon * ks .* besselj(1, epsilon * d * ks);
    M = MassMatrix(d, alpha, epsilon, L, N);
    Fs = M \ (omegas .* bs + epsilon^2 * gs / (alpha * sqrt(L)));
    Q = (sqrt(L) * d + 2 * dot(as, Gammas_d)) / (2 * sqrt(L) - 2 * dot(Fs, Gammas));
    if Q < 0
        Q
        warning("Q is negative");
    end
    
    %% Saves derivatives
    dydd(1 : N) = Q * Fs;
    dydd(N + 1 : 2 * N) = - Q * as;
end


end