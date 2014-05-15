function timePlotQ(a_N, a_Ystop)
global nsteps nnodes Xvec

% fID = fopen('Qdata1.dat','r');
% Qdata1 = fread(fID,[nsteps, nnodes]);
% fclose(fID);
% 
% fID = fopen('Qdata2.dat','r');
% Qdata2 = fread(fID,[nsteps,nnodes]);
% fclose(fID);

Qdata1 = dlmread('Qdata1.dat');
Qdata2 = dlmread('Qdata2.dat');

Q = cat(1,Qdata1,Qdata2);

for i=1:a_N*nsteps
    h = plot(Xvec,Q(i,:));
    ylim([-a_Ystop a_Ystop]);
    refreshdata(h);
    drawnow;
end

end