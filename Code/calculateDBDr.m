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

P = getPressure(a_A,a_R0,a_P0);
A0= pi*a_R0*a_R0;

%numTerm1 = (2.0*pi*a_R0/rho)*(P - a_P0);
%denTerm1 = 1.0 - (3.0*a_R0/(4.0*Elast*wallH))*(P - a_P0);

%P
%a_P0

%u = (A0/rho)*(P-a_P0);
%v = sqrt(A0/a_A);

% rho
% Elast
% wallH

%dudr = (((2*pi*a_R0)/rho)*(P-a_P0))-((4*Elast*wallH*pi)/(3*rho));
%dvdr = sqrt(pi/a_A);

%if abs(denTerm1) < 10^-30
%pause;
%end


%numTerm2 = (3.0*A0/(4.0*Elast*wallH*rho))*(P - a_P0)*(P - a_P0);
%denTerm2 = denTerm1^2;
%if abs(denTerm2) < 10^-30
%    pause;
%end
if ( isMaterialExpModel == 1 )
    f = (4.0/3.0)*(expModelK1*exp(expModelK2*a_R0) + expModelK3);
    dBdR0 = -f*pi/rho;
else
    dBdR0 = -(4.0*Elast*wallH*pi)/(3.0*rho);
end
%dBdR0 = ((dudr*v)-(u*dvdr))/(v*v)

%dBdR0 = (numTerm1/denTerm1) + (numTerm2/denTerm2);