function plotImpedance(omegas,zTree)
x0 = floor(length(omegas)/2)+1;

plot(omegas(x0:end),zTree(x0:end));
%axis([0 125 0 10^6]);
xlabel('Omega (2*pi*f)');
ylabel('Impedance (Z)');

end