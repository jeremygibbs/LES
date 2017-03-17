%%%% Created by Rob Stoll for ME EN 7960 008 LES of Turb Flow Spring 2009 
%%%% program to look at energy spectra and filtering of a 3D velocity field
%%%% currently using isotropic turbulence.  It uses the file iso_vel128.mat
%%%% containing the 3D U,V and W velocity fields (128^3 points).  NOTE:
%%%% current version assumes 2 filter sizes will be computed
clear all;close all;clc;

load iso_vel128.mat

DNS_bins = 128; %number of bins for Energy Spectra of DNS data
f = [32,16];    %number of points to filter over in filter subroutines (f=L/delta)

%%%% for the following options = 1 for yes (otherwise no)
opt1 = 1; %plot 3D spectra of original DNS data
opt2 = 0; %plot filtered 3D spectra
opt3 = 0; %plot 2D filtered velocity slices
opt4 = 0; %filter in physical space
opt5 = 0; %filter in spectral space
optp = 1; %output tiff file of plots

N=size(U);

if(opt1 == 1 || opt2 == 1)
    [E,DNS_BNS]=spectra3D(U,V,W,DNS_bins); %find the 3D Energy spectra
end
if(opt1 == 1)
    figure(1);loglog(DNS_BNS,E,'-ok','Linewidth',2,'Markersize',6);
    hold on;X1=5;X2=35;H=line([X1 X2],3*[X1^(-5/3) X2^(-5/3)]);
    hold on;L1=line([DNS_BNS(3) DNS_BNS(3)], get(gca,'YLim'));
    hold on;L2=line([DNS_BNS(45) DNS_BNS(45)], get(gca,'YLim'));
    text(15,0.05,'-5/3','FontSize',18);
    t2=text(15,0.05,'production','FontSize',18);
    t3=text(15,0.05,'isotropic scaling','FontSize',18);
    t4=text(15,0.05,'dissipation','FontSize',18);
    set(H,'color',[0 0 0],'LineWidth',2);
    set(L1,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    set(L2,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    set(t4, 'rotation', 90);
    axis([0.8 55 0.0001 0.4]);set(gca,'Fontsize',18)
    xlabel('k=(k_1^2+k_2^2+k_3^2)^{1/2}','Fontsize',18);ylabel('E(k)','Fontsize',18);
    if(optp == 1)
        print(gcf,'-dtiffnocompression','DNS_spectra');
    end
end

%%%% filtering the velocity field
if (opt5==1)
    display('calculating spectral-cutoff filter')
    [U_hat1,U_hat2] = cutoff_filt_3D(U,f);
    [V_hat1,V_hat2] = cutoff_filt_3D(V,f);
    [W_hat1,W_hat2] = cutoff_filt_3D(W,f);
    
    display('calculating gaussian filter')
    [U_ghat1,U_ghat2] = gaussian_filt_3D(U,f);
    [V_ghat1,V_ghat2] = gaussian_filt_3D(V,f);
    [W_ghat1,W_ghat2] = gaussian_filt_3D(W,f);
    
    display('calculating box filter')
    [U_bhat1,U_bhat2] = box_filt_3D(U,f);
    [V_bhat1,V_bhat2] = box_filt_3D(V,f);
    [W_bhat1,W_bhat2] = box_filt_3D(W,f);
end

if (opt4==1)
    display('calculating box filter')
	[U_bhat3,U_bhat4] = box_filt_3D_real(U,f);
	[V_bhat3,V_bhat4] = box_filt_3D_real(V,f);
	[W_bhat3,W_bhat4] = box_filt_3D_real(W,f);
    
    display('calculating Gaussian filter')
    [U_bhat5,U_bhat6] = gaussian_filt_3D_real(U,f);
    [V_bhat5,V_bhat6] = gaussian_filt_3D_real(V,f);
    [W_bhat5,W_bhat6] = gaussian_filt_3D_real(W,f);
    X=[0:N(1)-1]*2*pi/N(1);
    figure(7);pcolor(X,X,squeeze(U(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('DNS velocity','Fontsize',18)
    figure(8);pcolor(X,X,squeeze(U_bhat3(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('cutoff filter at \Delta_1','Fontsize',18)
    figure(9);pcolor(X,X,squeeze(U_bhat4(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('cutoff filter at \Delta_2','Fontsize',18)
    if(optp == 1)
        print(7,'-dtiffnocompression','DNS_2D');
        print(8,'-dtiffnocompression','cutoff_2D_D1');
        print(9,'-dtiffnocompression','cutoff_2D_D2');
    end
end

if(opt2 == 1 && opt5==1)
    display('computing box filter spectra')
    [E1,DNS_BNS]=spectra3D(U_hat1,V_hat1,W_hat1,DNS_bins); %find the 3D Energy spectra (cutoff)
    [E2,DNS_BNS]=spectra3D(U_hat2,V_hat2,W_hat2,DNS_bins); %find the 3D Energy spectra (cutoff)

    display('computing box filter spectra')
    [E1g,DNS_BNS]=spectra3D(U_ghat1,V_ghat1,W_ghat1,DNS_bins); %find the 3D Energy spectra (gaussian)
    [E2g,DNS_BNS]=spectra3D(U_ghat2,V_ghat2,W_ghat2,DNS_bins); %find the 3D Energy spectra (gaussian)

    [E1b,DNS_BNS]=spectra3D(U_bhat1,V_bhat1,W_bhat1,DNS_bins); %find the 3D Energy spectra (box)
    [E2b,DNS_BNS]=spectra3D(U_bhat2,V_bhat2,W_bhat2,DNS_bins); %find the 3D Energy spectra (box)

    figure(2);loglog(DNS_BNS,E,'-ok','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E1,'-sg','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E2,'-^r','Linewidth',2,'Markersize',6);
    hold on;Hk1=line([f(1) f(1)],[1e-7 10]);
    hold on;Hk2=line([f(2) f(2)],[1e-7 10]);
    set(Hk1,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    set(Hk2,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    axis([0.8 55 0.00001 0.4]);set(gca,'Fontsize',18)
    xlabel('k=(k_1^2+k_2^2+k_3^2)^{1/2}','Fontsize',18);ylabel('E(k)','Fontsize',18);
    title('Spectral cutoff filter','Fontsize',18)
    legend('DNS','\Delta_1','\Delta_2','Location','SouthWest')
    if(optp == 1)
        print(gcf,'-dtiffnocompression','cutoff_spectra');
    end
    
    figure(3);loglog(DNS_BNS,E,'-ok','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E1g,'-sg','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E2g,'-^r','Linewidth',2,'Markersize',6);
    hold on;Hk1=line([f(1) f(1)],[1e-7 10]);
    hold on;Hk2=line([f(2) f(2)],[1e-7 10]);
    set(Hk1,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    set(Hk2,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    axis([0.8 55 0.00001 0.4]);set(gca,'Fontsize',18)
    xlabel('k=(k_1^2+k_2^2+k_3^2)^{1/2}','Fontsize',18);ylabel('E(k)','Fontsize',18);
    title('Gaussian filter','Fontsize',18)
    legend('DNS','\Delta_1','\Delta_2','Location','SouthWest')
    if(optp == 1)
        print(gcf,'-dtiffnocompression','gaussian_spectra');
    end
    
    figure(4);loglog(DNS_BNS,E,'-ok','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E1b,'-sg','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E2b,'-^r','Linewidth',2,'Markersize',6);
    hold on;Hk1=line([f(1) f(1)],[1e-7 10]);
    hold on;Hk2=line([f(2) f(2)],[1e-7 10]);
    set(Hk1,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    set(Hk2,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    axis([0.8 55 0.00001 0.4]);set(gca,'Fontsize',18)
    xlabel('k=(k_1^2+k_2^2+k_3^2)^{1/2}','Fontsize',18);ylabel('E(k)','Fontsize',18);
    title('Box filter','Fontsize',18)
    legend('DNS','\Delta_1','\Delta_2','Location','SouthWest')
    if(optp == 1)
        print(gcf,'-dtiffnocompression','box_spectra');
    end
    
    figure(5);loglog(DNS_BNS,E,'-xk','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E1,'-or','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E1g,'-sg','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E1b,'-^b','Linewidth',2,'Markersize',6);
    hold on;Hk1=line([f(1) f(1)],[1e-7 10]);
    set(Hk1,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    axis([0.8 55 0.00001 0.4]);set(gca,'Fontsize',18)
    xlabel('k=(k_1^2+k_2^2+k_3^2)^{1/2}','Fontsize',18);ylabel('E(k)','Fontsize',18);
    title('Filter comparison at \Delta_1','Fontsize',18)
    legend('DNS','cutoff','gaussian','box','Location','SouthWest')
    if(optp == 1)
        print(gcf,'-dtiffnocompression','compare_spectra_D1');
    end
    
    figure(6);loglog(DNS_BNS,E,'-xk','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E2,'-or','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E2g,'-sg','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E2b,'-^b','Linewidth',2,'Markersize',6);
    hold on;Hk1=line([f(2) f(2)],[1e-7 10]);
    set(Hk1,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    axis([0.8 55 0.00001 0.4]);set(gca,'Fontsize',18)
    xlabel('k=(k_1^2+k_2^2+k_3^2)^{1/2}','Fontsize',18);ylabel('E(k)','Fontsize',18);
    title('Filter comparison at \Delta_2','Fontsize',18)
    legend('DNS','cutoff','gaussian','box','Location','SouthWest')
    if(optp == 1)
        print(gcf,'-dtiffnocompression','compare_spectra_D2');
    end
end

if (opt2==1 && opt4==1)
    [E3b,DNS_BNS]=spectra3D(U_bhat3,V_bhat3,W_bhat3,DNS_bins); %find the 3D Energy spectra (box)
    [E4b,DNS_BNS]=spectra3D(U_bhat4,V_bhat4,W_bhat4,DNS_bins); %find the 3D Energy spectra (box)
    [E5b,DNS_BNS]=spectra3D(U_bhat5,V_bhat5,W_bhat5,DNS_bins); %find the 3D Energy spectra (box)
    [E6b,DNS_BNS]=spectra3D(U_bhat6,V_bhat6,W_bhat6,DNS_bins); %find the 3D Energy spectra (box)
    
    figure(7);loglog(DNS_BNS,E,'-ok','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E3b,'-sg','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E4b,'-^r','Linewidth',2,'Markersize',6);
    hold on;Hk1=line([f(1) f(1)],[1e-7 10]);
    hold on;Hk2=line([f(2) f(2)],[1e-7 10]);
    set(Hk1,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    set(Hk2,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    axis([0.8 55 0.00001 0.4]);set(gca,'Fontsize',18)
    xlabel('k=(k_1^2+k_2^2+k_3^2)^{1/2}','Fontsize',18);ylabel('E(k)','Fontsize',18);
    title('Box filter','Fontsize',18)
    legend('DNS','\Delta_1','\Delta_2','Location','SouthWest')
    if(optp == 1)
        print(gcf,'-dtiffnocompression','box_real');
    end
    
    figure(8);loglog(DNS_BNS,E,'-ok','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E5b,'-sg','Linewidth',2,'Markersize',6);
    hold on;loglog(DNS_BNS,E6b,'-^r','Linewidth',2,'Markersize',6);
    hold on;Hk1=line([f(1) f(1)],[1e-7 10]);
    hold on;Hk2=line([f(2) f(2)],[1e-7 10]);
    set(Hk1,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    set(Hk2,'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    axis([0.8 55 0.00001 0.4]);set(gca,'Fontsize',18)
    xlabel('k=(k_1^2+k_2^2+k_3^2)^{1/2}','Fontsize',18);ylabel('E(k)','Fontsize',18);
    title('Gaussian filter','Fontsize',18)
    legend('DNS','\Delta_1','\Delta_2','Location','SouthWest')
    if(optp == 1)
        print(gcf,'-dtiffnocompression','gaussian_real');
    end
end
    
if(opt3 == 1)
    X=[0:N(1)-1]*2*pi/N(1);
    
    figure(9);pcolor(X,X,squeeze(U(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('DNS velocity','Fontsize',18)
    figure(10);pcolor(X,X,squeeze(U_hat1(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('cutoff filter at \Delta_1','Fontsize',18)
    figure(11);pcolor(X,X,squeeze(U_hat2(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('cutoff filter at \Delta_2','Fontsize',18)
    if(optp == 1)
        print(7,'-dtiffnocompression','DNS_2D');
        print(8,'-dtiffnocompression','cutoff_2D_D1');
        print(9,'-dtiffnocompression','cutoff_2D_D2');
    end
    
    figure(12);pcolor(X,X,squeeze(U_ghat1(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('gaussian filter at \Delta_1','Fontsize',18)
    figure(13);pcolor(X,X,squeeze(U_ghat2(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('gaussian filter at \Delta_2','Fontsize',18)
    if(optp == 1)
        print(10,'-dtiffnocompression','gaussian_2D_D1');
        print(11,'-dtiffnocompression','gaussian_2D_D2');
    end
    
    figure(14);pcolor(X,X,squeeze(U_bhat3(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('box filter at \Delta_1','Fontsize',18)
    figure(15);pcolor(X,X,squeeze(U_bhat4(:,:,64)));h=colorbar('vert');shading('interp');
    set(gca,'Fontsize',18);set(get(h,'Title'),'String','\fontsize{14} m/s');
    title('box filter at \Delta_2','Fontsize',18)
    if(optp == 1)
        print(12,'-dtiffnocompression','box_2D_D1');
        print(13,'-dtiffnocompression','box_2D_D2');
    end
end