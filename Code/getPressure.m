% Function to evaluate the vessle pressure assuming linear elastic response
% - a_A:    Value of area at the chosen spatial node
% - a_R0:   Value of initial radius at the chosen spatial node
% - a_P0:   Value of reference pressure for simulation
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------
function P = getPressure(a_A,a_R0,a_P0)
global Elast wallH             % mechanical properties of vessel wall
global expModelK1 expModelK2 expModelK3
global isMaterialExpModel

A0      = pi*a_R0*a_R0;

if ( isMaterialExpModel == 1 )
    f = (4.0/3.0)*(expModelK1*exp(expModelK2*a_R0) + expModelK3);
else
    f   = 4.0*Elast*wallH/(3.0*a_R0);
end

P       = a_P0 + f*(1.0 - sqrt(A0/a_A));

%if P < 0
%    P
%    error('Negative Pressure Value');
%end

