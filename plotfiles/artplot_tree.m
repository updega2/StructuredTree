clear all

txt1  = 'Aorta (1,5,7,9,11,13,15,17,19)';
txt2  = 'Brachiocephalic Artery (2)';
txt3  = 'Right Subclavian and Brachial Arteries (3)';
txt4  = 'Right Carotid Artery (4)';
txt6  = 'Left Carotid Artery (6)';
txt8  = 'Left Subclavian and Brachial Arteries (8)';
txt10 = 'Celiac Axis (10)';
txt12 = 'Superior Mesenteric (12)';
txt14 = 'Right Renal Artery (14)';
txt16 = 'Left Renal Artery (16)';
txt18 = 'Inferior Mesenteric (18)';
txt20 = 'External Iliac and Femoral Arteries (20,21,23)';
txt22 = 'Internal Iliac Artery (22)';
txt24 = 'Deep Femoral Arteries (24)';

files = [1 2 3 4 6 8 10 12 14 16 18 20 22 24];

zmi(1)  = -50;   % 1
zma(1)  = 450;   % 1
zmi(2)  =  -5;   % 2
zma(2)  = 100;   % 2
zmi(3)  = -20;   % 3
zma(3)  =  70;   % 3
zmi(4)  =  -2;   % 4
zma(4)  =  50;   % 4
zmi(5)  =  -2;   % 6
zma(5)  =  50;   % 6
zmi(6)  = -15;   % 8
zma(6)  =  75;   % 8
zmi(7)  = -20;   %10
zma(7)  = 225;   %10
zmi(8)  =  -2;   %12
zma(8)  =  45;   %12
zmi(9)  =  -2;   %14
zma(9)  =  27;   %14
zmi(10) =  -2;   %16
zma(10) =  27;   %16
zmi(11) =  -0.5; %18 
zma(11) =   4.25; %18
zmi(12) = -15;   %20
zma(12) =  55;   %20
zmi(13) =  -0.5; %22
zma(13) =   6.5; %22
zmi(14) = -10;   %24
zma(14) =  35;   %24

N=1024;  

figure(1); clf;
load -ascii ../runs/q1.2d;
[X,Y,Z] = gnuplot(q1); 

 h = plot(X(:,1),Z(:,1));

set(gca,'FontSize',18);
legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 -50 450]);
title('Aorta','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_asc.eps';

figure(2); clf;
load -ascii ../runs/p1.2d;
[Xp,Yp,Zp] = gnuplot(p1); 
h = plot(Xp(:,1),Zp(:,1));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Aorta','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_asc.eps';

figure(21); clf;
Zpm = mean(Zp);
h = plot(Yp(1,:),Zpm,Xp(:,32),Zp(:,32),Xp(:,64),Zp(:,64),Xp(:,96),...
    Zp(:,96),Xp(:,128),Zp(:,128),Xp(:,160),Zp(:,160));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
legend('mean',' 8 cm','16 cm','24 cm','32 cm','40 cm');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Aorta','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_steap.eps';

figure(22); clf;
Zpm = mean(Zp);
h = plot(Yp(1,:),Zpm);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',22);
axis([0 42 102 104.5]);
grid on;
set(gca,'YTick',[102.0 102.5  103.0 103.5 104.0 104.5]);
xlabel('x [cm]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_mean.eps';

print -depsc2 'p_aorta_per.eps'

figure(3); clf;
eps = 25;

h = plot(Xp(:,29),Z(:,29)-eps);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 -75 350]);
title('Aortic Arc','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_arc.eps';

figure(4); clf;

h = plot(Xp(:,29),Zp(:,29));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Aortic Arc','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_arc.eps';

figure(5); clf;
eps1 = 0.01;
eps2 = 29;

h = plot(Xp(:,39)+eps1,Z(:,39)-eps2);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 -50 300]);
title('Aortic Des','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_des.eps';

figure(6); clf;

h = plot(Xp(:,39),Zp(:,39));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Aortic Des','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_des.eps';

figure(7); clf;
eps=3;

h = plot(Xp(:,168),Z(:,168)-eps);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 -30 120]);
title('Aorta Bif','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_bif.eps';

figure(8); clf;

h = plot(Xp(:,168),Zp(:,168));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Aorta Bif','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_bif.eps';

figure(9)
eps1=0.02;
eps2=5;
load -ascii ../runs/q2.2d;
[X,Y,Z] = gnuplot(q2);

h = plot(X(:,8)+eps1,Z(:,8)+eps2);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 -15 90]);
title('Anonyma','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_ano.eps';

figure(10)
load -ascii ../runs/p2.2d;
[Xp,Yp,Zp] = gnuplot(p2);

h = plot(Xp(:,8),Zp(:,8));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Anonyma','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_ano.eps';

figure(11)
load -ascii ../runs/q4.2d;
[X,Y,Z] = gnuplot(q4);

eps2=2.2;
eps1=0.005;

h = plot(X(:,8),Z(:,8)-eps2);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 0 35]);
title('L Carotid','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_car.eps';

figure(12)
load -ascii ../runs/p4.2d;
[Xp,Yp,Zp] = gnuplot(p4);

h = plot(Xp(:,8),Zp(:,8));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('L Carotid','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_car.eps';

figure(13)
load -ascii ../runs/q8.2d;
[X,Y,Z] = gnuplot(q8);

eps2 = 6;
eps1 = 0.005;

h = plot(X(:,5)-eps1,Z(:,5)+eps2);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 -5 30]);
title('L Subclavian','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_sub.eps';

figure(14)
load -ascii ../runs/p8.2d;
[Xp,Yp,Zp] = gnuplot(p8);

h = plot(Xp(:,5),Zp(:,5));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('L Subclavian','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');


figure(15)
eps1 = 2;
eps2 = 0.025;
load -ascii ../runs/q20.2d;
[X,Y,Z] = gnuplot(q20);

h =plot(X(:,1)+eps2,Z(:,1)-eps1);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 -20 65]);
title('Iliac','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_ili.eps';

figure(16)
load -ascii ../runs/p20.2d;
[Xp,Yp,Zp] = gnuplot(p20);

h =plot(Xp(:,1),Zp(:,1));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Iliac','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_ili.eps';

figure(17)
eps=2;

h = plot(Xp(:,32),Z(:,32)-eps);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 -20 50]);
title('Femoral','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_fem.eps';

figure(18)

h = plot(Xp(:,32),Zp(:,32));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Femoral','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_fem.eps';

figure(19)
eps = 0.9;
load -ascii ../runs/q3.2d;
[X,Y,Z] = gnuplot(q3);

h = plot(X(:,134),Z(:,134)-eps);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 -1 7]);
title('Brachials','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_bra.eps';

figure(20)
load -ascii ../runs/p3.2d;
[Xp,Yp,Zp] = gnuplot(p3);

h = plot(Xp(:,134),Zp(:,134));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);

legend('Simulations');
set(gca,'FontSize',22);
axis([3.25 4.3 80 140]);
title('Brachials','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_bra.eps';




