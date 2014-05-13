function kcfl = checkCFL(a_A,a_Q)
global delK delH
global nnodes
global Xvec

kcfl = 100;

for i=1:nnodes
    c = getWavespeed(Xvec(i),a_A(i));
    V = a_Q(i)/a_A(i);
    
    temp = min(delH/(abs(V-c)),delH/(abs(V+c)));
    kcfl = min(temp,kcfl);
end
if delK > kcfl
    error('CFL violated, minimum k required is %f',kcfl);
end
end

    