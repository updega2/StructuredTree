% Function to evaluate the derivative of source terms with respect to the
% independent flow variables.
% Argument definitions:
% - a_A:    Value of area at the chosen spatial node
% - a_Q:    Value of flow rate at the chosen spatial node
% - a_R0:   Value of initial radius at the chosen spatial node
% - a_dR0dX:Value of the derivative of the initial radius distribution
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------
function DSDU = jacobianSource(a_A,a_Q,a_R0,a_dR0dX)
%Define global variables
global Elast wallH wallDelta    % mechanical properties of vessel wall
global Nu rho                   % physical properties of the fluid


DSDU = ones(2,2)*10^-30;

A0  = pi*a_R0*a_R0;
C1  = 2.0*pi*a_R0/rho;
C2  = 3.0*a_R0/(4.0*Elast*wallH);
C3  = 3.0*A0/(4.0*Elast*wallH*rho);
R   = sqrt(a_A/pi);

PSI         = (4.0*Elast*wallH/(3.0*a_R0))*(1 - sqrt(A0/a_A));
dPSIdA      = -(4.0*Elast*wallH/(3.0*a_R0))*sqrt(A0)*(-0.5)*(a_A^(-3.0/2.0));
dBdR0_dA    = (C1/((1 - C2*PSI)^2) + 2.0*C3*PSI/((1.0 - C2*PSI)^3))*dPSIdA;
dBdR0_dA    = dBdR0_dA*a_dR0dX;

DSDU(1,1)   = 0.0;
DSDU(1,2)   = 0.0;
DSDU(2,1)   = (2.0*pi*Nu*R/wallDelta)*(a_Q/(a_A^2)) + dBdR0_dA;
DSDU(2,2)   = -2.0*pi*Nu*R/(wallDelta*a_A);



