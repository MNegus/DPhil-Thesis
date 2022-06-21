function [Es_outer, Es_splash_sheet] = dropletenergy(ts, TimeDependents, SubstrateCoefficients, epsilon)
%DROPLETENERGY Summary of this function goes here
%   Detailed explanation goes here
    
    ds = TimeDependents.ds;
    d_ts = TimeDependents.d_ts;

    Es_outer = outerenergy(ds, d_ts, SubstrateCoefficients, epsilon);
    Es_splash_sheet = splashsheetenergy(ts, ds, d_ts, epsilon);
end

