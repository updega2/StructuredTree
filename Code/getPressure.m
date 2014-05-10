% Function to evaluate the vessle pressure assuming linear elastic response
% - a_A:    Value of area at the chosen spatial node
% - a_R0:   Value of initial radius at the chosen spatial node
% - a_P0:   Value of reference pressure for simulation
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------
function P = getPressure(a_A,a_R0,a_P0)
global Elast wallH             % mechanical properties of vessel wall

A0      = pi*a_R0*a_R0;
coeff   = 4.0*Elast*wallH/(3.0*a_R0);

P       = a_P0 + coeff*(1.0 - sqrt(A0/a_A));
