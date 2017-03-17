function [tau,Coeff] = Burgers_SGS(u,dudx,dx,opt)
%lam = 2 for explicit filtering. use lam = 1 for delta/h = 1

[R C]   = size(u); if R<C u = u'; end
N       = length(u);
u2      = [u(N);u;u(1)];
u4      = [u(N-1);u(N);u;u(1);u(2)];
Coeff   = 0;
%*******************************************************************************
if  opt == 1 
    lam     = 1;                                                %Constant Coeff Smagorinsky
    CS2     = 0.16^2;
    
    d1      = dealias1(abs(dudx),N);
    d2      = dealias1(dudx,N);
    d3      = dealias2(d1.*d2,N);
    
    tau     = -2*CS2*((lam*dx)^2)*d3;
    Coeff   = sqrt(CS2);
%**************************************************************************
%*******************************************************************************
elseif  opt == 2 
   % Place your 1st model here!
%**************************************************************************
elseif  opt == 3 
   % Place your 2nd model here!
%**************************************************************************
elseif  opt == 4 
   % Place your 3rd model here!
%**************************************************************************
else                                                            %NO SGS Model
    tau     = zeros(size(u)); Coeff = 0;
end

%*******************************************************************************    
function x_p= dealias1(x,N)
    M       = N/2;
    fx      = fft(x);
    fx_p    = [fx(1:N/2+1)' zeros(1,M) fx(N/2+2:N)']';
    x_p     = real(ifft(fx_p));
return
%*******************************************************************************    
function x   = dealias2(x_p,N)
    M        = N/2;
    fx_p     = fft(x_p);
    fx       = [fx_p(1:N/2+1)' fx_p(M+N/2+2:M+N)']';
    fx(N/2+1)= 0;
    x        = (3/2)*real(ifft(fx));
return
%Last Updated: 24th September, 2005***************************************************