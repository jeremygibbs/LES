function x1 = Burgers_FBM(alpha,m)
%This program generates the noise term for the Stochastic Burgers equation
%Written by: Sukanta Basu************************************************* 
x           = sqrt(m)*randn(1,m);
L           = m/2 ;
freq        = [ 1 1:L (L-1):-1:1 ];
xft         = fft(x);
xft(1)      = 0; xft(L+1) = 0;  %Set the Nyquist to zero
x1ft        = xft .* (freq .^ (- alpha/2 ) );
x1          = real( ifft(x1ft) );
%Last Updated: November 13th, 2003****************************************