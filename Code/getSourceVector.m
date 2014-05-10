function Svec = getSourceVector(a_A,a_Q,a_R0,a_dR0dX,a_P0)
%--------------------------------------------------------------------------
% Function to obtain the source vector for the hyperbolic set of flow
% equations presented in the paper.
% Argument definition
% - a_A:    Value of area at the chosen spatial node
% - a_R0:   Value of initial radius at the chosen spatial node
% - a_P0:   Value of reference pressure for simulation
% - a_Q:    Value of flow rate at the chosen spatial node
% - a_dR0dX:Value of the derivative of the initial radius distribution
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------

Svec = zeros(2,1);

if (a_A <0)
    a_A
    error('A_a is negative');
end

Svec(1) = 0.0;
Svec(2) = wallShear1(a_A,a_Q) + calculateDBDr(a_A,a_R0,a_P0)*a_dR0dX;