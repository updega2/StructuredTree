function [] = plotFlowRateIn(a_N, a_CardOut, a_CardMax, a_CardPeriod)
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

t0 = 0.0;
t1 = a_N*a_CardPeriod;
tVec = linspace(t0,t1);
qVec = ones(size(tVec))*10^-30;

for k1 = 1:length(tVec)
    qVec(k1) = getFlowRateIn(tVec(k1),a_CardOut, a_CardMax, a_CardPeriod);
end

figure()
plot(tVec, qVec,'r','linewidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','bold');
ylabel('Cardiac Flow Rate (cc/s)','FontSize',12,'FontWeight','bold');