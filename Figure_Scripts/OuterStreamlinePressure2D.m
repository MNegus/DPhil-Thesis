%% OuterStreamlinePressure
% Script to produce a figure of the streamlines and pressure in the outer
% region for the 2D impact case. Appears in Section 3.3.3, under the Wagner
% theory chapter. 
% 
% For the figure, we neglect any motion of the membrane, taking w == 0. In
% this case, the complex potential is given by
% W(\zeta) =  \phi + i \psi = i(\zeta^2 - d(t)^2)^{1/2}, 
% and the pressure is
% p = \Re(i A(t) / (\zeta^2 - d(t)^2)^{1/2}), A(t) = d(t) d'(t).
%
% We plot the streamlines as a contour plot for constant \psi, and then a
% colour plot for the pressure p. We pick a generic value of d(t) for
% plotting clarity. 