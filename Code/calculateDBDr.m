function dBdR0 = calculateDBDr(a_A,a_R0,a_P0)
%--------------------------------------------------------------------------
% Function to evaluate the derivative of B with respect to R0, as defined
% in the paper by Olufsen and Peskin
% Argument definitions:
% - a_A:    Value of area at the chosen spatial node
% - a_R0:   Value of initial radius at the chosen spatial node
% - a_P0:   Value of reference pressure for simulation
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------
global Elast wallH           % mechanical properties of vessel wall
global rho                   % physical property of the fluid
global isMaterialExpModel    % choice variable for Olufsen's exponential material model
global expModelK1 expModelK2 expModelK3

P   = getPressure(a_A,a_R0,a_P0);
A0  = pi*a_R0*a_R0;
R   = sqrt(a_A/pi);


if ( isMaterialExpModel == 1 )
    f       = (4.0/3.0)*(expModelK1*exp(expModelK2*a_R0) + expModelK3);
    dfdR0   = (4.0/3.0)*(expModelK1*expModelK2)*exp(expModelK2*a_R0);
else
    f       = (4.0*Elast*wallH)/(3.0*a_R0);
    dfdR0   = -(4.0*Elast*wallH)/(3.0*(a_R0*a_R0));
end

term1 = dfdR0*((pi*a_R0*R-pi*a_R0*a_R0)/rho);
term2 = (f/rho)*(pi*R-2*pi*a_R0);
dBdR0 = term1+term2;
