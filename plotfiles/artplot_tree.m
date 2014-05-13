clear all

load met06sor_ano01.DAT
load met06sor_arc01.DAT
load met06sor_asc01.DAT
load met06sor_bif01.DAT
load met06sor_bra01.DAT
load met06sor_car01.DAT
load met06sor_des01.DAT
load met06sor_fem01.DAT
load met06sor_ili01.DAT
load met06sor_sub01.DAT

t_ano = met06sor_ano01(:,1)./1000+3*1.10;
q_ano = met06sor_ano01(:,2);
t_arc = met06sor_arc01(:,1)./1000+3*1.10;
q_arc = met06sor_arc01(:,2);
t_asc = met06sor_asc01(:,1)./1000+3*1.10;
q_asc = met06sor_asc01(:,2);
t_bif = met06sor_bif01(:,1)./1000+3*1.10;
q_bif = met06sor_bif01(:,2);
t_bra = met06sor_bra01(:,1)./1000+3*1.10;
q_bra = met06sor_bra01(:,2);
t_car = met06sor_car01(:,1)./1000+3*1.10;
q_car = met06sor_car01(:,2);
t_des = met06sor_des01(:,1)./1000+3*1.10;
q_des = met06sor_des01(:,2);
t_fem = met06sor_fem01(:,1)./1000+3*1.10;
q_fem = met06sor_fem01(:,2);
t_ili = met06sor_ili01(:,1)./1000+3*1.10;
q_ili = met06sor_ili01(:,2);
t_sub = met06sor_sub01(:,1)./1000+3*1.10;
q_sub = met06sor_sub01(:,2);

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

%for k=1:length(files)
%  zmi(k) = 80;
%  zma(k) = 160;
%end;
%%zmi(13) = 70;
%%zma(13) = 180;

N=1024;  

figure(1); clf;
load -ascii q1.2d;
[X,Y,Z] = gnuplot(q1); 
%load -ascii q1.2d;
%[Xnd,Ynd,Znd] = gnuplot(q1); 
x = t_asc;
y = q_asc;
xi=(0:N-1)'/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
size(x)
size(y)
size(Z(:,1))
size(xi)
%h = plot(xi,yi,xi,Z(:,1),xi,Znd(:,1));
h = plot(xi,yi,xi,Z(:,1));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim ');
set(gca,'FontSize',22);
axis([3.25 4.3 -50 450]);
title('Aorta','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_asc.eps';

figure(2); clf;
load -ascii p1.2d;
[Xp,Yp,Zp] = gnuplot(p1); 
%load -ascii p1.2d;
%[Xp1,Yp1,Zp1] = gnuplot(p1); 
%h = plot(xi,Zp(:,1),xi,Zp1(:,1));
h = plot(xi,Zp(:,1));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Aorta','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_asc.eps';

figure(21); clf;
Zpm = mean(Zp);
h = plot(Yp(1,:),Zpm,xi,Zp(:,32),xi,Zp(:,64),xi,Zp(:,96),xi,Zp(:,128),xi,Zp(:,160));
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

%figure(23); clf;
%load -ascii p1_long.2d;
%[Xp,Yp,Zp] = gnuplot(p1_long); 
%h = plot(Xp(:,32),Zp(:,32));
%for i=1:length(h)
% get(h(i));
% set(h(i),'LineWidth',3);
%end;
%set(gca,'FontSize',22);
%axis([0 8 50 130]);
%grid on;
%xlabel('t [s]');
%ylabel('p [mmHg]');
%print -depsc2 'p_aorta_per.eps'

figure(3); clf;
eps = 25;
x = t_arc;
y = q_arc;
xi=(0:N-1)'/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
%h = plot(xi,yi,xi,Z(:,29)-eps,xi,Znd(:,29)-eps);
h = plot(xi,yi,xi,Z(:,29)-eps);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim');
set(gca,'FontSize',22);
axis([3.25 4.3 -75 350]);
title('Aortic Arc','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_arc.eps';

figure(4); clf;
%h = plot(xi,Zp(:,29),xi,Zp1(:,29));
h = plot(xi,Zp(:,29));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
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
x = t_des;
y = q_des;
xi=(0:N-1)'/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
%h = plot(xi,yi,xi+eps1,Z(:,39)-eps2,xi+eps1,Znd(:,39)-eps2);
h = plot(xi,yi,xi+eps1,Z(:,39)-eps2);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim');
set(gca,'FontSize',22);
axis([3.25 4.3 -50 300]);
title('Aortic Des','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_des.eps';

figure(6); clf;
%h = plot(xi,Zp(:,39),xi,Zp1(:,39));
h = plot(xi,Zp(:,39));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Aortic Des','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_des.eps';

figure(7); clf;
eps=3;
x = t_bif;
y = q_bif;
xi=(0:N-1)'/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
%h = plot(xi,yi,xi,Z(:,168)-eps,xi,Znd(:,168)-eps); % FIT WK
%h = plot(xi,yi,xi,Z(:,168)-eps,xi,Znd(:,168)-5*eps); % FIT WK
h = plot(xi,yi,xi,Z(:,168)-eps);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim');
set(gca,'FontSize',22);
axis([3.25 4.3 -30 120]);
title('Aorta Bif','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_bif.eps';

figure(8); clf;
%h = plot(xi,Zp(:,168),xi,Zp1(:,168));
h = plot(xi,Zp(:,168));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
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
load -ascii q2.2d;
[X,Y,Z] = gnuplot(q2);
%load -ascii q2.2d;
%[Xnd,Ynd,Znd] = gnuplot(q2);
x = t_ano;
y = q_ano;
xi=(0:N-1)/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
%h = plot(xi,yi,xi+eps1,Z(:,8)+eps2,xi+eps1,Znd(:,8)+eps2);
h = plot(xi,yi,xi+eps1,Z(:,8)+eps2);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim ');
set(gca,'FontSize',22);
axis([3.25 4.3 -15 90]);
title('Anonyma','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_ano.eps';

figure(10)
load -ascii p2.2d;
[Xp,Yp,Zp] = gnuplot(p2);
%load -ascii p2.2d;
%[Xp1,Yp1,Zp1] = gnuplot(p2);
%h = plot(xi,Zp(:,8),xi,Zp1(:,8));
h = plot(xi,Zp(:,8));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Anonyma','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_ano.eps';

figure(11)
load -ascii q4.2d;
[X,Y,Z] = gnuplot(q4);
%load -ascii q4.2d;
%[Xnd,Ynd,Znd] = gnuplot(q4);
eps2=2.2;
eps1=0.005;
x = t_car;
y = q_car;
xi=(0:N-1)/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
%h = plot(xi,yi,xi-eps1,Z(:,8)-eps2,xi-eps1,Znd(:,8)-eps2);
h = plot(xi,yi,xi-eps1,Z(:,8)-eps2);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim');
set(gca,'FontSize',22);
axis([3.25 4.3 0 35]);
title('L Carotid','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_car.eps';

figure(12)
load -ascii p4.2d;
[Xp,Yp,Zp] = gnuplot(p4);
%load -ascii p4.2d;
%[Xp1,Yp1,Zp1] = gnuplot(p4);
%h = plot(xi,Zp(:,8),xi,Zp1(:,8));
h = plot(xi,Zp(:,8));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('L Carotid','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_car.eps';

figure(13)
load -ascii q8.2d;
[X,Y,Z] = gnuplot(q8);
%load -ascii q8.2d;
%[Xnd,Ynd,Znd] = gnuplot(q8);
eps2 = 6;
eps1 = 0.005;
x = t_sub;
y = q_sub;
xi=(0:N-1)/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
%h = plot(xi,yi,xi-eps1,Z(:,5)+eps2,xi-eps1,Znd(:,5)+eps2);
h = plot(xi,yi,xi-eps1,Z(:,5)+eps2);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim');
set(gca,'FontSize',22);
axis([3.25 4.3 -5 30]);
title('L Subclavian','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_sub.eps';

figure(14)
load -ascii p8.2d;
[Xp,Yp,Zp] = gnuplot(p8);
%load -ascii p8.2d;
%[Xp1,Yp1,Zp1] = gnuplot(p8);
%h = plot(xi,Zp(:,5),xi,Zp(:,5));
h = plot(xi,Zp(:,5));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('L Subclavian','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_sub.eps';

figure(15)
eps1 = 2;
eps2 = 0.025;
load -ascii q20.2d;
[X,Y,Z] = gnuplot(q20);
%load -ascii q20.2d;
%[Xnd,Ynd,Znd] = gnuplot(q20);
x = t_ili;
y = q_ili;
xi=(0:N-1)/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
%h =plot(xi,yi,xi+eps2,Z(:,1)-eps1,xi+eps2,Znd(:,1)-4*eps1);  % FIT WK
%h =plot(xi,yi,xi+eps2,Z(:,1)-eps1,xi+eps2,Znd(:,1)-4*eps1);
h =plot(xi,yi,xi+eps2,Z(:,1)-eps1);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim');
set(gca,'FontSize',22);
axis([3.25 4.3 -20 65]);
title('Iliac','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_ili.eps';

figure(16)
load -ascii p20.2d;
[Xp,Yp,Zp] = gnuplot(p20);
%load -ascii p20.2d;
%[Xp1,Yp1,Zp1] = gnuplot(p20);
%h =plot(xi,Zp(:,1),xi,Zp1(:,1));
h =plot(xi,Zp(:,1));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Iliac','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_ili.eps';

figure(17)
eps=2;
x = t_fem;
y = q_fem;
xi=(0:N-1)/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
%h = plot(xi,yi,xi,Z(:,32)-eps,xi,Znd(:,32)-2*eps);  % FIT WK
%h = plot(xi,yi,xi,Z(:,32)-eps,xi,Znd(:,32)-eps);
h = plot(xi,yi,xi,Z(:,32)-eps);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim');
set(gca,'FontSize',22);
axis([3.25 4.3 -20 50]);
title('Femoral','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_fem.eps';

figure(18)
%h = plot(xi,Zp(:,32),xi,Zp1(:,32));
h = plot(xi,Zp(:,32));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
set(gca,'FontSize',22);
axis([3.25 4.3 80 130]);
title('Femoral','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_fem.eps';

figure(19)
eps = 0.9;
load -ascii q3.2d;
[X,Y,Z] = gnuplot(q3);
%load -ascii q3.2d;
%[Xnd,Ynd,Znd] = gnuplot(q3);
x = t_bra;
y = q_bra;
xi=(0:N-1)/N*(x(end)-x(1)) + x(1); 
yi = spline(x,y,xi);
%h = plot(xi,yi,xi,Z(:,134)-eps,xi,Znd(:,134)-eps);
h = plot(xi,yi,xi,Z(:,134)-eps);
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('meas','sim st','sim wk');
legend('meas','sim');
set(gca,'FontSize',22);
axis([3.25 4.3 -1 7]);
title('Brachialis','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('q [cm^3/s]');
print -depsc2 'aorta_bra.eps';

figure(20)
load -ascii p3.2d;
[Xp,Yp,Zp] = gnuplot(p3);
%load -ascii p3.2d;
%[Xp1,Yp1,Zp1] = gnuplot(p3);
%h = plot(xi,Zp(:,134),xi,Zp1(:,134));
h = plot(xi,Zp(:,134));
for i=1:length(h)
 get(h(i));
 set(h(i),'LineWidth',3);
end;
set(gca,'FontSize',18);
%legend('sim st','sim wk');
legend('sim');
set(gca,'FontSize',22);
axis([3.25 4.3 80 140]);
title('Brachialis','Fontsize',24);
grid on;
xlabel('t [s]');
ylabel('p [mmHg]');
print -depsc2 'p_aorta_bra.eps';




