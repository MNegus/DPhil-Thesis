function [M, S] = MassMatrix(d, alpha, epsilon, L, N)
%MASSMATRIX Outputs the mass matrix for the normal modes solution

    %% Save lambdas
    ks = pi * (2 * (1 : N) - 1) / (2 * L);

    %% Calculate the S matrix
    % Bulk nodes for m ~= n
    S = (pi * d / epsilon) ...
        * (besselj(0, epsilon * d * ks)' * (ks .* besselj(1, epsilon * d * ks)) ...
        - (ks .* besselj(1, epsilon * d * ks))' * besselj(0, epsilon * d * ks)) ...
        ./ (ks.^2 - (ks.^2)');
    
    % Diagonals for m == n
    S(1 : N + 1 : end) = (pi * d^2 / 2) ...
        * (besselj(0, epsilon * d * ks).^2 + besselj(1, epsilon * d * ks).^2);
    
    %% Write the overall matrix
    M = eye(N) + epsilon^2 * S / (alpha * L);
    
end
