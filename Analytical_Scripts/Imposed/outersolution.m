function [Ws, ps] = outersolution(zetas, t, SubstrateFunctions)
%OUTERSOLUTION Returns the leading-order solution in the outer region
%   Returns the leading-order solution for W, the complex potential, and p,
%   the pressure, in the outer region for a quadratic substrate.

    %% Only works for the 2D case
    if SubstrateFunctions.dimension ~= '2D'
       error("Invalid dimension. Only dimension == '2D' is supported.");
    end

    %% Load substrate coefficients
    d = SubstrateFunctions.d(t);
    A = SubstrateFunctions.A(t);
    
    a_t = SubstrateFunctions.a_t(t);
    a_tt = SubstrateFunctions.a_tt(t);
    
    b_t = SubstrateFunctions.b_t(t);
    b_tt = SubstrateFunctions.b_tt(t);
    
    %% Complex potential
    % Complex integral solutions
    I1_0 = 1 - zetas ./ (zetas.^2 - d^2).^(1/2);
    I1_1 = d^2 / 2 + zetas.^2 - zetas.^3 ./ (zetas.^2 - d^2).^(1/2);
    
    % Solution for Ws, the complex potential
    Ws = 1i * (zetas.^2 - d^2).^(1/2) .* (1 - a_t * I1_0 - (b_t / 3) * I1_1);
    
    %% Pressure
    % Complex integral solutions
    I3_0 = zetas .* ((zetas.^2 - d^2).^(1/2) - zetas) + d^2 / 2;
    I3_1 = (1 / 8) * (4 * d^2 * zetas.^2 ...
        + 8 * zetas.^3 .* ((zetas.^2 - d^2).^(1/2) - zetas) + d^4);
    
    % Solution for ps, the pressure
    ps = real(1i * (A - a_tt * I3_0 - (b_tt / 3) * I3_1) ./ (zetas.^2 - d^2).^(1/2));
end

