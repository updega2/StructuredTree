function [] = plotPressure(a_N, a_x_loc,a_R0, a_P0,a_CardPeriod)
%--------------------------------------------------------------------------
% Function to visualize the inflow flow rate for the one dimensional blood
% flow model as presented by Olufsen and Peskin
% Argument Definition:
% a_N:          The number of cardiac cycles to be plotted
% a_CardOut:    Total cardiac output (in cc/sec)
% a_CardMax:    The instant of maximal cardiac output (in sec)
% a_CardPeriod: The period of cardiac cycle (in sec)
% Authors: Debanjan Mukherjee, Adam Updegrove, Alex Baelde
% University of California, Berkeley.
%--------------------------------------------------------------------------
global solA

t0 = 0.0;
t1 = a_N*a_CardPeriod;
tVec = linspace(t0,t1);
pVec = ones(size(tVec))*10^-30;

for k1 = 1:length(tVec)
    pVec(k1) = getPressure(solA(k1,a_x_loc),a_R0(a_x_loc), a_P0);
end

figure()
plot(tVec, pVec,'r','linewidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','bold');
ylabel('Pressure (Pa)','FontSize',12,'FontWeight','bold');