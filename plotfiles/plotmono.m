clear all

N=1024;  

figure(1); clf;
load -ascii q1.2d;
[X,Y,Z] = gnuplot(q1); 
h = plot(X(:,39),Z(:,39));
set(h,'LineWidth',3);
set(gca,'FontSize',18);
%axis([3.25 4.3 -50 450]);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_des.eps'

figure(2); clf;
load -ascii p1.2d;
[Xp,Yp,Zp] = gnuplot(p1); 
h = plot(Xp(:,39),Zp(:,39));
set(h,'LineWidth',3);
set(gca,'FontSize',18);
%axis([3.25 4.3 80 130]);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_asc.eps'
