function hs = jetfreesurface(xs, t, SubstrateFunctions)

    %% Find w
    if SubstrateFunctions.dimension == "2D"
        error("Currently only supported for dimension == 'axi'.");
    else
        w = SubstrateFunctions.w(t);
    end

    %% Load functions
    epsilon = SubstrateFunctions.epsilon;
    d = SubstrateFunctions.d;
    d_t = SubstrateFunctions.d_t;
    d_tt = SubstrateFunctions.d_tt;
    J = SubstrateFunctions.J;
    
    %% Find minimum tau (from manual xMax)
    xMax = max(xs);
    xBarMax = xMax / epsilon;
    zeroFun = @(tau) xBarMax - 2 * d_t(tau) * (t - tau) - d(tau);
    tauMin = fsolve(zeroFun, 1e-6);
    
    %% Find minumum tau (from setting minimum hBar)
%     hBarMin = 1e-4;
%     zeroFun = @(tau) hBarMin - (d_t(tau) .* J(tau)) ./ (d_t(tau) - 2 * d_tt(tau) .* (t - tau));
%     tauMin = fsolve(zeroFun, 1e-6);

    %% Find taus, clustered near t
    lambda = 10;
    b = (tauMin - t) / (1 - exp(-lambda));
    a = tauMin - b;
    taus = a + b * exp(-lambda * linspace(0, 1, 1e3));
    
%     taus = linspace(tauMin, t, 1e3);

    %% Determinant of Jacobian, throw error if zero at any point
    detJ = d_t(taus) - 2 * d_tt(taus) .* (t - taus);
    if (sum(detJ <= 0) > 1)
        detJ
        error("Singular Jacobian");
    end

    %% Find jet free surface
    ds = SubstrateFunctions.d(taus);
    d_ts = SubstrateFunctions.d_t(taus);
    d_tts = SubstrateFunctions.d_tt(taus);
    Js = SubstrateFunctions.J(taus);

    xBars = 2 * d_ts .* (t - taus) + ds;
    hBars = (d_ts .* Js) ./ (d_ts - 2 * d_tts .* (t - taus));
    uBars = 2 * d_ts; % Velocity, if needed
    
    %% Determine solution in original variables
%     xs = epsilon * xBars;
    hs = -epsilon^2 * w + epsilon^3 * interp1(epsilon * xBars, hBars, xs, ...
        'linear', 'extrap');

end