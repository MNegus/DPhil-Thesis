function [pMaxs, xMaxs] = pressuremax(ts, SubstrateFunctions)
%PRESSUREMAX Summary of this function goes here
%   Detailed explanation goes here

    epsilon = SubstrateFunctions.epsilon;
    ds = SubstrateFunctions.d(ts);
    d_ts = SubstrateFunctions.d_t(ts);
    Js = SubstrateFunctions.J(ts);
    
    pMaxs = d_ts.^2 / (2 * epsilon^2);
    xMaxs = epsilon * ds - epsilon^3 * 6 * Js / pi; 
end

