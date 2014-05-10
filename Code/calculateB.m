function B = calculateB(a_A,a_R0,a_P0)
%--------------------------------------------------------------------------
% Function to evaluate the term B as defined in paper by Olufsen and Peskin
% Argument definitions:
% - a_A:    Value of area at the chosen spatial node
% - a_R0:   Value of initial radius at the chosen spatial node
% - a_P0:   Value of reference pressure for simulation
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------
global Elast wallH              % mechanical properties of vessel wall
global rho                    % physical property of the fluid
global isMaterialExpModel
global expModelK1 expModelK2 expModelK3

%P   = getPressure(a_A,a_R0,a_P0);
A0  = pi*a_R0*a_R0;

%numTerm = (A0/rho)*(P - a_P0);
%denTerm = 1.0 - (3.0*a_R0/(4.0*Elast*wallH))*(P - a_P0);

if ( isMaterialExpModel == 1 ) 
    f = (4.0/3.0)*(expModelK1*exp(expModelK2*a_R0) + expModelK3);
else
    f = (4.0*Elast*wallH)/(3.0*a_R0);
end

% if abs(denTerm) < 10^-30
%     pause;
% end

B = ((sqrt(A0*a_A)*f)/rho)*(1-sqrt(A0/a_A));
%B = numTerm/denTerm;




