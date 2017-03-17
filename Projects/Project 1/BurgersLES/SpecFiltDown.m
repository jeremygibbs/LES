function [uf] = SpecFiltDown(u,k)
%Fourier Filtering (k-Delta)
%Written by: Sukanta Basu**************************************************
N       = length(u);
M       = N/k;
fu      = fft(u);

fuf             = zeros(M,1);
fuf(1:M/2)      = fu(1:M/2);
fuf(M/2+2:M)    = fu(N-M/2+2:N);

uf              = (1/k)*real(ifft(fuf));
