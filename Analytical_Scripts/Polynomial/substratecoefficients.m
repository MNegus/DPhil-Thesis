function SubstrateCoefficients ...
    = substratecoefficients(as, bs, a_ts, b_ts, a_tts, b_tts, epsilon)
%%substratecoefficients
% Loads the coefficients of the quadratic substrate, such that
% w(x, t) = a(t) + b(t) * x^2,
% where we define aHat = a, bHat = epsilon^2 * b for the outer solution
% part.
    
    % a coefficients
    SubstrateCoefficients.aHats = as;
    SubstrateCoefficients.aHat_ts = a_ts;
    SubstrateCoefficients.aHat_tts = a_tts;
    
    % b coefficients
    SubstrateCoefficients.bHats = epsilon^2 * bs;
    SubstrateCoefficients.bHat_ts = epsilon^2 * b_ts;
    SubstrateCoefficients.bHat_tts = epsilon^2 * b_tts;
end