function timePlotA(t,n,ystop)
global solA

for i=1:t
    h = plot(linspace(1,n,n),solA(i,:));
    ylim([0 ystop]);
    refreshdata(h);
    drawnow;
end

end