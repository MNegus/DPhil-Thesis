function N = NMax(alpha, beta, gamma, L, q, epsilon, delta_t)
%NMAX Finds the maximum allowable value of N for a given delta_t

    if (gamma > 0)
        c = alpha * gamma * (2 * pi / (q * delta_t))^2;
        lambda2 = (sqrt(epsilon^4 * beta^2 + 4 * c) - epsilon^2 * beta) / (2 * gamma);
        N = 0.5 + (L / (epsilon * pi)) * sqrt(lambda2);
    else
        error("gamma = 0 not supported.")
        N = 0.5 + (2 * L / pi) * sqrt(alpha / beta);
    end
end