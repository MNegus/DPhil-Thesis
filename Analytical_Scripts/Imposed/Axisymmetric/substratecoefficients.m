function SubstrateCoefficients ...
    = substratecoefficients(ws, w_ts, w_tts)
%%substratecoefficients
% Loads the coefficients of the quadratic substrate, such that
% w(x, t) = w(t).
    
    SubstrateCoefficients.ws = ws;
    SubstrateCoefficients.w_ts = w_ts;
    SubstrateCoefficients.w_tts = w_tts;
end