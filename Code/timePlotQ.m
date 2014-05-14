function timePlotQ(t,n,ystop)
global solQ

for i=1:t
    h = plot(linspace(1,n,n),solQ(i,:));
    ylim([0 ystop]);
    refreshdata(h);
    drawnow;
end

end