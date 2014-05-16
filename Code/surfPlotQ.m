function surfPlotQ(a_N, a_Ystop)
global nsteps nnodes Xvec P0 Tcard

Qdata1 = dlmread('Qdata1.dat');
Qdata2 = dlmread('Qdata2.dat');
Qdata3 = dlmread('Qdata3.dat');
Qdata4 = dlmread('Qdata4.dat');

Qdata1 = Qdata1.*10^6;
Qdata2 = Qdata2.*10^6;
Qdata3 = Qdata3.*10^6;
Qdata4 = Qdata4.*10^6;

Q = cat(1,Qdata1,Qdata2,Qdata3,Qdata4);
Tvec = linspace(0, Tcard, nsteps);
Tvec = Tvec + Tcard*ones(size(Tvec));
Tvec = Tvec + Tcard*ones(size(Tvec));
Tvec = Tvec + Tcard*ones(size(Tvec));

figure()
X = surf(Xvec,Tvec,Qdata1,'LineStyle','none');
xlabel('x [cm]');
ylabel('time [sec]');
zlabel('q [cm^3/s]');
set(gca, 'YDir', 'reverse')
cmap = contrast(X);
colormap(cmap);

figure()
Tvec = Tvec + Tcard*ones(size(Tvec));
Y = surf(Xvec,Tvec,Qdata2,'Linestyle','none');
xlabel('x [cm]');
ylabel('time [sec]');
zlabel('q [cm^3/s]');
set(gca, 'YDir', 'reverse')
cmap = contrast(Y);
colormap(cmap);

figure()
Tvec = Tvec + Tcard*ones(size(Tvec));
Z = surf(Xvec,Tvec,Qdata3,'LineStyle','none');
xlabel('x [cm]');
ylabel('time [sec]');
zlabel('q [cm^3/s]');
set(gca, 'YDir', 'reverse')
cmap = contrast(Z);
colormap(cmap);

figure()
Tvec = Tvec + Tcard*ones(size(Tvec));
K = surf(Xvec,Tvec,Qdata4,'LineStyle','none');
xlabel('x [cm]');
ylabel('time [sec]');
zlabel('q [cm^3/s]');
set(gca, 'YDir', 'reverse')
cmap = contrast(K);
colormap(cmap);

end