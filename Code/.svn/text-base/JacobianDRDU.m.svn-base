% Function to evaluate the Jacobian Matrix of the flux terms with respect
% to the driving variables.
% Argument definitions:
% - a_A:    Value of area at the chosen spatial node
% - a_R0:   Value of initial radius at the chosen spatial node
% - a_P0:   Value of reference pressure for simulation
% - a_Q:    Value of flow rate at the chosen spatial node
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------
function DRDU = JacobianDRDU(a_A,a_Q,a_R0,a_P0)
%Define global variables
global Elast wallH              % mechanical properties of vessel wall
global rho                      % physical property of the fluid

DRDU = ones(2,2)*10^-30;

P   = getPressure(a_R0,a_A,a_P0);
A0  = pi*a_R0*a_R0;

Dpmp0Da  = -(4.0*Elast*wallH/(3.0*a_R0))*sqrt(A0)*(-0.5)*(a_A^(-3.0/2.0));
numTerm1 = (A0/rho)*Dpmp0Da*(1.0 - ...
    (3.0*a_R0/(4.0*Elast*wallH))*(P - a_P0));
numTerm2 = (A0/rho)*(P - a_P0)*(-3.0*a_R0/(4.0*Elast*wallH))*Dpmp0Da;
denTerm  = 1.0 - (3.0*a_R0/(4.0*Elast*wallH))*(P - a_P0);
dBdA     = (numTerm1 - numTerm2)/denTerm;
DRDU(1,1) = 0.0;
DRDU(1,2) = 1.0;
DRDU(2,1) = -(a_Q/a_A)*(a_Q/a_A) + dBdA;
DRDU(2,2) = 2.0*a_Q/a_A;