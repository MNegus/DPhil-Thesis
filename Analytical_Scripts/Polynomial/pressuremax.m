function [pMaxs, xMaxs] = pressuremax(TimeDependents, epsilon)
%PRESSUREMAX Summary of this function goes here
%   Detailed explanation goes here


    ds = TimeDependents.ds;
    d_ts = TimeDependents.d_ts;
    Js = TimeDependents.Js;
    
    pMaxs = d_ts.^2 / (2 * epsilon^2);
    xMaxs = epsilon * ds - epsilon^3 * 6 * Js / pi; 
end

