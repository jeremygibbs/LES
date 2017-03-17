%This Program Solves Stochastic Burgers Equation Using*********************** 
%the Fourier Collocation Method**********************************************
clear all; clc; close all; t = cputime;

%Input Parameters************************************************************
%Nx:    # Grid Points; dx: Grid Spacing; NSteps: Total # Iterations
%dt:    Time Increment; visc = Kinematic Viscosity; D: Diffusivity of Noise
Nx  = 8192; dx = 2*pi/Nx; NSteps = 2e6; dt = 1e-4; 
visc = 1e-5; D = 1e-6; 

outpath = ['.',filesep,'BurgersDNS_',int2str(Nx),filesep];
[~,~,~] = mkdir(outpath);
file_u  = fopen([outpath,'DNS_',int2str(Nx),'.out'],'w');

%Initial velocity field******************************************************
randn('state',0);                   %Initialize Random Number Generator
u           = zeros(Nx,1); fu = fft(u); fu(Nx/2+1) = 0; u = real(ifft(fu));
%Advancing in time***********************************************************
for s = 1:NSteps
    %Compute the derivatives*************************************************
    [dudx,d2udx2,du2dx]  = Burgers_Derivatives(u,dx,1);
    %Create Noise ***********************************************************
    f   = Burgers_FBM(0.75,Nx)';
    
    %Right Hand Side*********************************************************
    RHS     = visc*d2udx2 - 0.5*du2dx + sqrt(2*D/dt)*f;
    if s == 1
        unew = u + dt*RHS;                              %Euler
    else    
        unew = u + dt*(1.5*RHS - 0.5*RHSp);             %AB2
    end
    funew = fft(unew); funew(Nx/2+1) = 0; unew = real(ifft(funew));
    u     = unew; RHSp = RHS;
    
    if sum(isnan(u))>0 break; end;
    
    if rem(s,100) == 0
        CFL         = max(abs(u))*dt/dx;            
        KE          = 0.5*var(u);
        t           = cputime-t; 
        fprintf('%d\t%f\t%f\t%f\t%f\t%f\t%f\n',s,s*dt,KE,CFL,max(u),min(u),t);
        t           = cputime;
    end
    
    %Output Space-Time Data**************************************************       
    if rem(s,1000) == 0
        fprintf(file_u,'%f\t',u); fprintf(file_u,'\n');
    end
end%End of s loop
fclose(file_u); 
%Last Updated: September 24, 2005**************************************************