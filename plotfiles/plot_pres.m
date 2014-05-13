clear all

N=1024; 

figure(1);clf;
%load -ascii ../runs/p6.2d;
load -ascii ../runs/p4.2d;
[Xp4,Yp4,Zp4] = gnuplot(p4); % Should be changed to 6 for right side
%load -ascii ../runs/out1/p6.2d;
%load -ascii ../runs/out1/p4.2d; % Should be changed to 6 for right side
load -ascii ../runs/p4.2d
[Xp4a,Yp4a,Zp4a] = gnuplot(p4);
L4 = size(Zp4,2);
h = plot(Xp4(:,1),Zp4(:,L4),Xp4(:,1),Zp4a(:,L4));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
legend('constr','normal');
set(gca,'FontSize',22);
axis([3.25 4.5 75 130]);
title('Pressure in carotids');
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_lc.eps'

figure(2)
%load -ascii ../runs/p8.2d;
load -ascii ../runs/p3.2d; 
[Xp3,Yp3,Zp3] = gnuplot(p3); % Should be changed to 8 for right side
%load -ascii ../runs/out1/p8.2d;
%load -ascii ../runs/out1/p3.2d;
load -ascii ../runs/p3.2d;
[Xp3a,Yp3a,Zp3a] = gnuplot(p3); % Should be changed to 8 for right side
L3 = size(Zp3,2);
h = plot(Xp3(:,1),Zp3(:,L3),Xp3(:,1),Zp3a(:,L3));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
legend('constr','normal');
set(gca,'FontSize',22);
axis([3.25 4.5 75 130]);
title('Pressure in brachials');
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_lb.eps'

figure(3)
h = plot(Xp3(:,1),Zp3(:,L3),Xp3(:,1),Zp4(:,L4));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
legend('carotids','brach');
set(gca,'FontSize',22);
axis([3.25 4.5 75 130]);
title('Comparison of Pressures');
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_lbc.eps'

figure(4)
%load -ascii ../runs/q6.2d;
load -ascii ../runs/q4.2d;
[Xq1,Yq1,Zq1] = gnuplot(q4); % Should be changed to 6 for right side
%load -ascii ../runs/out1/q6.2d;
%load -ascii ../runs/out1/q4.2d;
load -ascii ../runs/q4.2d;
[Xq,Yq,Zq] = gnuplot(q4); % Should be changed to 6 for right side
h = plot(Xq1(:,1),Zq1(:,L4),Xq1(:,1),Zq(:,L4));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
legend('constr','normal');
set(gca,'FontSize',22);
axis([3.25 4.5 4 20]);
title('Flow in carotids');
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'q_lc.eps'

figure(5)
%load -ascii ../runs/q8.2d;
load -ascii ../runs/q3.2d;
[Xq1,Yq1,Zq1] = gnuplot(q3); % Should be changed to 8 for right side
%load -ascii ../runs/out1/q8.2d;
%load -ascii ../runs/out1/q3.2d;
load -ascii ../runs/q3.2d;
[Xq,Yq,Zq] = gnuplot(q3); % Should be changed to 8 for right side
h = plot(Xq1(:,1),Zq1(:,L3),Xq1(:,1),Zq(:,L3));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
legend('constr','normal');
set(gca,'FontSize',22);
axis([3.25 4.5 0.75 4]);
title('Flow in brachials');
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'q_lb.eps'







