function Y = oFFTshift(Z)

N = int64(length(Z));

Y = circshift(Z,N/int64(2));

end