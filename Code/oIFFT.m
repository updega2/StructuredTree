function Z = oIFFT(Zhat)

N = max(size(Zhat));

M = round(log(real(N))/log(2.));

Z = Zhat;

for l = 1:M
    le = 2^l;
    le1 = le/2;
    u =complex(1,0);
    w = complex(cos(pi,le1),sin(pi/le1));
    
    for j = 1:le1
        for i = j:le:N
            ip = 1+le1;
            t = Z(ip)*u;
            Z(i) = Z(i) + t;
        end
        u = u*w;
    end
end

end

            