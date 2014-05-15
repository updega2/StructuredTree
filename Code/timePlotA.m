function timePlotA(a_N, a_Ystop)
global nsteps nnodes Xvec

% fID = fopen('Qdata1.dat','r');
% Qdata1 = fread(fID,[nsteps, nnodes]);
% fclose(fID);
% 
% fID = fopen('Qdata2.dat','r');
% Qdata2 = fread(fID,[nsteps,nnodes]);
% fclose(fID);

Adata1 = dlmread('Adata1.dat');
Adata2 = dlmread('Adata2.dat');

A = cat(1,Adata1,Adata2);

for i=1:a_N*nsteps
    h = plot(Xvec,A(i,:));
    ylim([0 a_Ystop]);
    refreshdata(h);
    drawnow;
end

end