function Q1 = getFlowRateIn(a_T, a_CardOut, a_CardMax, a_CardPeriod)
%--------------------------------------------------------------------------
% Function to calculate the inflow flow rate for the one dimensional blood
% flow model as presented by Olufsen and Peskin
% Argument Definition:
% a_T:          The time at which the flow rate is sought
% a_CardOut:    Total cardiac output (in m^3/sec)
% a_CardMax:    The instant of maximal cardiac output (in sec)
% a_CardPeriod: The period of cardiac cycle (in sec)
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------

t   = mod(a_T, a_CardPeriod);
Q1  = ((a_CardOut*t)/(a_CardMax^2))*exp(-(t^2)/(2*(a_CardMax)^2));

