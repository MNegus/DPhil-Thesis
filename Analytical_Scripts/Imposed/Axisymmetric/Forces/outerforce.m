function Fs = outerforce(ds, d_ts, d_tts, epsilon)
%%outerforce
%   Returns the outer solution for the force on the substrate in the
%   axisymmetric, plate impact case.

    Fs = (8 * epsilon / 9) * ds.^3 .* (ds .* d_tts + 4 * d_ts.^2);
    
    % At t == 0, set the force to be 0
    zeroIdx = find(ds == 0);
    if ~isempty(zeroIdx)
       Fs(zeroIdx) = 0;
    end
end