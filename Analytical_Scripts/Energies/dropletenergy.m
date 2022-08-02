function [Es_outer, Es_jets] = dropletenergy(ts, SubstrateFunctions)
%DROPLETENERGY Summary of this function goes here
%   Detailed explanation goes here

    Es_outer = outerenergy(ts, SubstrateFunctions);
    Es_jets = jetsenergy(ts, SubstrateFunctions);
end

