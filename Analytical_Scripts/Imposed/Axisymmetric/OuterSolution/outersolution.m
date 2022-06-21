function [Ws, ps] = outersolution(zetas, d, A, SubstrateCoefficients)
%OUTERSOLUTION Returns the leading-order solution in the outer region
%   Returns the leading-order solution for W, the complex potential, and p,
%   the pressure, in the outer region for a quadratic substrate.

    %% Load substrate coefficients
    aHat_t = SubstrateCoefficients.aHat_ts;
    aHat_tt = SubstrateCoefficients.aHat_tts;
    
    bHat_t = SubstrateCoefficients.bHat_ts;
    bHat_tt = SubstrateCoefficients.bHat_tts;
    
    %% Complex potential
    % Complex integral solutions
    I1_0 = 1 - zetas ./ (zetas.^2 - d^2).^(1/2);
    I1_1 = d^2 / 2 + zetas.^2 - zetas.^3 ./ (zetas.^2 - d^2).^(1/2);
    
    % Solution for Ws, the complex potential
    Ws = 1i * (zetas.^2 - d^2).^(1/2) .* (1 - aHat_t * I1_0 - (bHat_t / 3) * I1_1);
    
    %% Pressure
    % Complex integral solutions
    I3_0 = zetas .* ((zetas.^2 - d^2).^(1/2) - zetas) + d^2 / 2;
    I3_1 = (1 / 8) * (4 * d^2 * zetas.^2 ...
        + 8 * zetas.^3 .* ((zetas.^2 - d^2).^(1/2) - zetas) + d^4);
    
    % Solution for ps, the pressure
    ps = real(1i * (A - aHat_tt * I3_0 - (bHat_tt / 3) * I3_1) ./ (zetas.^2 - d^2).^(1/2));
end

