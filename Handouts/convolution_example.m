%%%% File does an example convolution using a spatial box filter in 1D %%%%
clear all;close all;clc;

f=[16]; %number of points to filter over (f=L/delta)

load iso_vel128.mat;

Ux = squeeze(U(:,64,64));
clear U;clear V;clear W;
xdim = length(Ux);
xf   = zeros(xdim);
xbig = [Ux(xdim-f/2:xdim);Ux;Ux(1:f/2)];          %periodic in x
xl = linspace(1,128,128);
xf=Ux;
v = VideoWriter('newfile.mp4','MPEG-4');
v.FrameRate = 8;
v.Quality = 100;
open(v);
for i=1:xdim
    xbox=zeros(xdim,1);
    xf(i)=mean(xbig(i:i+f));
    if((i-f) < 0)
        xbox(1:i)=1;
    elseif(i > xdim-f/2)
        xbox(i:xdim)=1;
    else
        xbox(i-f/2:i+f/2)=1;
    end
    figure(1);
    subplot(211);plot(xl,Ux,xl,xf,'LineWidth',2);axis([1 128 -1.5 1.5]);
    ylabel('$\phi$','Interpreter','LaTex','FontSize',24)
    ax = gca;ax.XTick=[0,32,64,96,128];
    subplot(212);plot(xbox,'LineWidth',2);axis([1 128 0 2]);
    xlabel('x','FontSize',24);ylabel('$G$','Interpreter','LaTex','FontSize',24)
    ax = gca;ax.XTick=[0,32,64,96,128];
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.1)  
end
close(v)