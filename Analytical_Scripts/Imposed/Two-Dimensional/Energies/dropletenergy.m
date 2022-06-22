function [Es_outer, Es_jets] = dropletenergy(ts, SubstrateFunctions, epsilon)
%DROPLETENERGY Summary of this function goes here
%   Detailed explanation goes here

    Es_outer = outerenergy(ts, SubstrateFunctions, epsilon);
    Es_jets = jetsenergy(ts, SubstrateFunctions, epsilon);
end

