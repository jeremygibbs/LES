function [dudx,d2udx2,du2dx,d3udx3] = Burgers_Derivatives(u,dx,opt)
%This program computes different order derivatives based on**************** 
%Fourier Collocation & Central Difference scheme of second or fourth order* 
%accuracy******************************************************************
%Written by: Sukanta Basu**************************************************
%Make sure the series is N rows and 1 column*******************************
[NR,NC] = size(u); if NC>1 u = u'; end 
N       = length(u);

%Fourier Collocation Method************************************************
if opt == 1
h       = 2*pi/N;
fac     = h/dx;
k       = [0 1:(N/2-1) 0 -(N/2-1):1:-1]';
fu      = fft(u);
dudx    = (fac)*real(ifft(sqrt(-1)*k.*fu));
d2udx2  = (fac^2)*real(ifft(-k.*k.*fu));
d3udx3  = (fac^3)*real(ifft(-sqrt(-1)*k.^3.*fu));
%d4udx4  = (fac^4)*real(ifft((k.^4).*fu));
%To compute du2dx we need dealiasing***************************************
%Zero padding for dealiasing***********************************************
fu_p    = [fu(1:N/2+1)' zeros(1,N) fu(N/2+2:N)']';
u_p     = real(ifft(fu_p));
u2_p    = u_p.*u_p;
fu2_p   = fft(u2_p);
fu2     = [fu2_p(1:N/2+1)' fu2_p(N+N/2+2:N+N)']';
du2dx   = (2)*(fac)*real(ifft(sqrt(-1)*k.*fu2));
%The factor (2) comes from {(2)^2}*(1/2) [due to changing from N to 2N]
end

%2nd order accuracy********************************************************
if opt == 2 
%We apply periodic boundary condition for padding**************************
up      = [u(N-1);u(N);u;u(1);u(2)]; 
vp      = up.^2;
i       = 3:1:length(up)-2;
%Derivatives***************************************************************
dudx    = (up(i+1) - up(i-1))/(2*dx); 
d2udx2  = (up(i+1) - 2*up(i) + up(i-1))/(dx^2);
du2dx   = (vp(i+1) - vp(i-1))/(2*dx);
%d4udx4  = (up(i+2) - 4*up(i+1) +6*up(i) - 4*up(i-1) + up(i-2))/(dx^4);
end

%4th order accuracy********************************************************
if opt == 4
%We apply periodic boundary condition for padding**************************
up      = [u(N-2);u(N-1);u(N);u;u(1);u(2);u(3)]; 
vp      = up.^2;
i       = 4:1:length(up)-3;
%Derivatives***************************************************************
dudx    = (- up(i+2) + 8*up(i+1) - 8*up(i-1) + up(i-2))/(12*dx); 
d2udx2  = (- up(i+2) + 16*up(i+1) - 30*up(i) + 16*up(i-1) - up(i-2))/(12*(dx^2));
du2dx   = (- vp(i+2) + 8*vp(i+1) - 8*vp(i-1) + vp(i-2))/(12*dx);
%d4udx4  = (- up(i+3) + 12*up(i+2) - 39*up(i+1) + 56*up(i) -39*up(i-1) + 12*up(i-2)...
%          - up(i-3))/(6*(dx^4));
end
%Last Updated: 12th Nov, 2003**********************************************
