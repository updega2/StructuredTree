function Rvec = getFluxVector(a_A,a_Q,a_R0,a_P0)
%--------------------------------------------------------------------------
% Function to obtain the flux vector for the hyperbolic set of flow
% equations defined in the paper.
% Argument definitions:
% - a_A:    Value of area at the chosen spatial node
% - a_R0:   Value of initial radius at the chosen spatial node
% - a_P0:   Value of reference pressure for simulation
% - a_Q:    Value of flow rate at the chosen spatial node
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------

Rvec = zeros(2,1);

if (a_A <0)
    a_A
    error('A_a is negative');
end

Rvec(1) = a_Q;
Rvec(2) = a_Q*a_Q/a_A + calculateB(a_A,a_R0,a_P0);

