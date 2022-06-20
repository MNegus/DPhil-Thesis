function [Es_outer, Es_jets] = dropletenergy(ts, TimeDependents, SubstrateCoefficients, epsilon)
%DROPLETENERGY Summary of this function goes here
%   Detailed explanation goes here

    Es_outer = outerenergy(TimeDependents.ds, SubstrateCoefficients, epsilon);
    Es_jets = jetsenergy(ts, TimeDependents.Bs, TimeDependents.Cs, epsilon);
end

