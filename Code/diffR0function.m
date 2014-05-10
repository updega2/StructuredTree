% Function to obtain the equilibrium vessel radius spatial derivative
% using the assumption of exponential taper as presented by 
% Olufsen and Peskin
% Argument definition
% - a_X:    Coordinate at the chosen spatial node
% - a_Rt:   Value of inlet (top) radius for the chosen vessel segment
% - a_Rb:   Value of outlet (bottom) radius for the chosen vessel segment
% - a_L:    Length of the chosen vessel segment
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------
function dR0dX = diffR0function(a_X, a_Rt, a_Rb, a_L)

R0X = getR0function(a_X,a_Rt,a_Rb,a_L);
dR0dX = R0X*((log(a_Rb/a_Rt))*(1.0/a_L));
%dR0dX = a_Rt*exp((log(a_Rb/a_Rt))*(a_X/a_L))*(1/a_L)*(log(a_Rb/a_Rt));
