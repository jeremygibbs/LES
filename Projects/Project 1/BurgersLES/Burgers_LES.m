%This Program Solves Stochastic Burgers Equation Using*********************** 
%the Fourier Collocation Method**********************************************
clear all; clc; close all; t = cputime;

%Input Parameters************************************************************
%Nx:    # Grid Points; dx: Grid Spacing; NSteps: Total # Iterations
%dt:    Time Increment; visc = Kinematic Viscosity; D: Diffusivity of Noise
Nx  = 512; dx = 2*pi/Nx; NSteps = 2e6; dt = 1e-4; 
visc = 1e-5; D = 1e-6; 

%A few files for writing output
outpath = ['.',filesep,'BurgersLES_',int2str(Nx),filesep];
[~,~,~] = mkdir(outpath);

file_u  = fopen([outpath,'Smag_',int2str(Nx),'_U.out'],'w');
file_Cs = fopen([outpath,'Smag_',int2str(Nx),'_Coeff.out'],'w');
fw      = fopen([outpath,'Smag_',int2str(Nx),'_SGS.out'],'w');

%Initial velocity field******************************************************
randn('state',0);                   %Initialize Random Number Generator
u           = zeros(Nx,1); fu = fft(u); fu(Nx/2+1) = 0; u = real(ifft(fu));
%Advancing in time***********************************************************
for s = 1:NSteps
    %Compute the derivatives*************************************************
    [dudx,d2udx2,du2dx,d3udx3]  = Burgers_Derivatives(u,dx,1);
    %Create Noise ***********************************************************
    f   = Burgers_FBM(0.75,8192)'; ff = SpecFiltDown(f,8192/Nx); %Note change FBM res (8192) to match DNS
    
    %Compute SGS Stress and its Divergence***********************************
    [tau,Coeff] = Burgers_SGS(u,dudx,dx,1);
    dtaudx      = Burgers_Derivatives(tau,dx,1); 
    
    %Right Hand Side*********************************************************
    RHS     = visc*d2udx2 - 0.5*du2dx + sqrt(2*D/dt)*ff -0.5*dtaudx;
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
        fprintf('%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',s,s*dt,KE,CFL,max(u),min(u),Coeff,t);
        t           = cputime;
    end
    
    %Output Space-Time Data**************************************************       
    if rem(s,1000) == 0
        fprintf(file_u,'%f\t',u); fprintf(file_u,'\n');
        fprintf(file_Cs,'%f\n',Coeff); 
        
        E1 = mean(-tau.*dudx);
        E2 = mean(visc*dudx.^2);
        E3 = mean(-tau.*d3udx3); 
        E4 = mean(dudx.^3);
        E5 = mean(visc*d2udx2.^2);
        fprintf(fw,'%d\t%f\t%f\t%f\t%f\t%f\n',i,E1,E2,E3,E4,E5);    
    end
    
end%End of s loop
fclose all;
%Last Updated: September 24, 2005**************************************************