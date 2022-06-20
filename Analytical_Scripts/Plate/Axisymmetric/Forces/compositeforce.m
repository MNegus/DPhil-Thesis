function Fs = compositeforce(Fs_outer, Fs_inner, Fs_overlap)
%%compositeforce
%   Returns the additive composite expansion of the force on the substrate
    Fs = Fs_outer + Fs_inner - Fs_overlap;
end