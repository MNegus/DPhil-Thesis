function Fs = outerforce(ds, As, w_tts)
%%outerforce
%   Returns the outer solution for the force on the substrate in the
%   two-dimensional, plate impact case.

    Fs = pi * As - w_tts .* ds.^3 / 3;
end