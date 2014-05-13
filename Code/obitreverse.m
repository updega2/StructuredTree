function Y = obitreverse(X)
int64 j;
int64 k;

N = int64(max(size(X)));
j=1;

for l=1:N-1
    
    Y(j) = X(l);
    k = idivide(N,int64(2));
    
    while k < j
        j = j-k;
        k = idivide(k,int64(2));
    end
    j= j+k;
    l;
end

Y(N) = X(N);

end
    