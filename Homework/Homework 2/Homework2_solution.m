%%% homework #2 mfile
close all;clear all;clc;warning off

%%% problem #1

%%% Read in velocity data
vel = load('iso_vel128.mat');
u   = vel.U;
v   = vel.V;
w   = vel.W;

%%% Get size of velocity data
nx  = size(u,1);
ny  = size(u,2);
nz  = size(u,3);

%%% Compute wavenumbers
kx  = 2*pi*((-0.5*nx):0.5*(nx))/nx; 
ky  = 2*pi*((-0.5*ny):0.5*(ny))/ny; 
kz  = 2*pi*((-0.5*nz):0.5*(ny))/nz;
k   = 2*pi*(0:0.5*nx)/nx;

kx  = (-0.5*nx):0.5*(nx); 
ky  = (-0.5*ny):0.5*(ny); 
kz  = (-0.5*nz):0.5*(ny);
mk  = ( (0.5*nx)^2 + (0.5*ny)^2 + (0.5*nz)^2)^0.5
k   = (1:mk)
%%% Compute spectra
uf  = fftshift(fftn(u)/nx/ny/nz);
vf  = fftshift(fftn(v)/nx/ny/nz);
wf  = fftshift(fftn(w)/nx/ny/nz);

uk = uf.*conj(uf);
vk = vf.*conj(vf);
wk = wf.*conj(wf);

%filter in real space
delta  = 1/16
ndelta = nx*delta

u_flt=zeros(nx,ny,nz);
v_flt=zeros(nx,ny,nz);
w_flt=zeros(nx,ny,nz);
wght=ones(ndelta+1,ndelta+1,ndelta+1);
for ii =1: nx
    for jj =1: ny
        for kk =1: nz
        normwght=wght(max(1,0.5*ndelta-ii+2):min(ndelta+1,nx-ii+1+0.5*ndelta), max(1,0.5*ndelta-jj+2):min(ndelta+1,ny-jj+1+0.5*ndelta),max(1,0.5*ndelta-kk+2):min(ndelta+1,nz-kk+1+0.5*ndelta));
        normwght=normwght/sum(sum(sum(normwght))); 
        u_flt(ii,jj,kk)=sum(sum(sum(normwght.*u(max(1,ii-0.5*ndelta):min(nx,ii+0.5*ndelta),max(1,jj-0.5*ndelta):min(nx,jj+0.5*ndelta),max(1,kk-0.5*ndelta):min(nx,kk+0.5*ndelta)))));
        v_flt(ii,jj,kk)=sum(sum(sum(normwght.*v(max(1,ii-0.5*ndelta):min(nx,ii+0.5*ndelta),max(1,jj-0.5*ndelta):min(nx,jj+0.5*ndelta),max(1,kk-0.5*ndelta):min(nx,kk+0.5*ndelta)))));
        w_flt(ii,jj,kk)=sum(sum(sum(normwght.*w(max(1,ii-0.5*ndelta):min(nx,ii+0.5*ndelta),max(1,jj-0.5*ndelta):min(nx,jj+0.5*ndelta),max(1,kk-0.5*ndelta):min(nx,kk+0.5*ndelta)))));
        end
    end
end

uf_flt  = fftshift(fftn(u_flt)/nx/ny/nz);
vf_flt  = fftshift(fftn(v_flt)/nx/ny/nz);
wf_flt  = fftshift(fftn(w_flt)/nx/ny/nz);
uk_flt = uf_flt.*conj(uf_flt);
vk_flt = vf_flt.*conj(vf_flt);
wk_flt = wf_flt.*conj(wf_flt);

Ev  = uk+vk+wk;
Ev_flt  = uk_flt+vk_flt+wk_flt;

Eb  = zeros(numel(k),1);
Eb_flt  = zeros(numel(k),1);

%%% Find energy in 
for ii =1: nx
    for jj =1: ny
        for kk =1: nz 
            kh=(kx(ii)^2+ky(jj)^2+kz(kk)^2)^0.5; 
            ind=find(k-kh>=0);
            if(~isempty(ind))
                Eb(ind(1))=Eb(ind(1))+Ev(ii,jj,kk);
                Eb_flt(ind(1))=Eb_flt(ind(1))+Ev_flt(ii,jj,kk);
            end
        end
    end
end
%%% Plot spectra
[m,ind] = max(Eb);
k_kol = 1.5*(m/k(ind)^-(5/3))*k.^-(5/3); 
k_isr = [k(4) k(41)];
f     = figure('Position',[100 100 600 450]); 
loglog(k,Eb,'k-','LineWidth',2)
hold on
loglog(k,Eb_flt,'k-','LineWidth',2)
hold on
loglog(k,k_kol,'b--','LineWidth',2) 
l1=line([k_isr(1) k_isr(1)],get(gca,'Ylim')); 
set(l1,'Color','r','LineStyle','--','LineWidth',1.5) 
l2=line([k_isr(2) k_isr(2)],get(gca,'Ylim')); 
set(l2,'Color','r','LineStyle','--','LineWidth',1.5)
grid on
leg=legend('$E(k)$','$E\propto k^{-5/3}$'); 
set(leg,'Interpreter','Latex','Location','SouthWest')
%xlim([0.95*k(1) 1.05*k(end)]) 
%ylim([1E-5 1E0]) 
xlabel('$k$','Interpreter','LaTex','FontSize',18) 
ylabel('$E(k)$','Interpreter','LaTex','FontSize',18)
print('3dspectra','-dpng')