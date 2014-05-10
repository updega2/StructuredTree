function C = getWavespeed(a_X, a_A)

% Include all the relevant global variables in the function at the
% beginning:
%--------------------------------------------------------------------------
global Elast wallH              % mechanical properties of vessel wall
global rho                      % physical property of the fluid
global isMaterialExpModel
global expModelK1 expModelK2 expModelK3
global isR0function
global Rtop Rbottom
global X1 X0

LenVessel = X1 - X0;

if ( isR0function )
    R0    = getR0function(a_X, Rtop, Rbottom, LenVessel);
else
    R0    = getR0data(Xm); %%%% NEED TO IMPLEMENT THIS FUNCTION
end

A0 = pi*R0*R0;

if ( isMaterialExpModel == 1 )
    f = (4.0/3.0)*(expModelK1*exp(expModelK2*a_R0) + expModelK3);
else
    f = (4.0*Elast*wallH)/(3.0*a_R0);
end

C = sqrt((0.5*f/rho)*sqrt(A0/a_A));