function tauRHS = wallShear1(a_A,a_Q)
%--------------------------------------------------------------------------
% Function to calculate the source terms due to wall shear stress as
% presented in the paper by Olufsen and Peskin
% - a_A:    Value of area at the chosen spatial node
% - a_Q:    Value of flow rate at the chosen spatial node
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------
global Nu wallDelta 
 
if (a_A <0)
    a_A
    error('A_a is negative');
end

R      = sqrt(a_A/pi);
tauRHS = - (2.0*pi*Nu*R*a_Q)/(wallDelta*a_A);