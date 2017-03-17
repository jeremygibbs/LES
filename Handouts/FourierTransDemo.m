%%% Example problem from Stull 88 (Ex. 8.4.2) %%%
%%% example wavelet (Haar) transform          %%%
close all;clear all;clc
ScreenSize=get(0,'Screensize'); %check the size of the screen we are plotting on
width=floor(ScreenSize(3)*4/5);height=floor(ScreenSize(4)*4/5); %dimensions of figure plot
%%% the data points %%%%%%
q=[8 9 9 6 10 3 5 6]; %space/time series data goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%
lnsp={'-k';'--b';':r';'-.g';'-c';'--m';':y';'-.k';'-b';'--r'}; %for plotting needs to be length(q)
N=length(q); %how many points
X=[0:500]*(N-1)/500; %array for plotting

%%% calculate the Fourier coefficients
%%% (note matlab doesn't normalize forward)

fk=fft(q)/N;

cosine=@(x,ak,n,N)ak*cos(2*pi*n*x/N);
sine=@(x,bk,n,N)bk*sin(2*pi*n*x/N);

wavesum=zeros(size(X));
figure(1);set(gcf,'Position',[floor(ScreenSize(3)*0.1) floor(ScreenSize(4)*0.55) width height])

figure(1);subplot(311);grid on;box on %real (cosine) part
title(['wave forms for k=']);ylabel('real(f_k)');
figure(1);subplot(312);grid on;box on %imaginary (sine) part
ylabel('imag(f_k)');
figure(1);subplot(313);hold off;grid on;box on
plot([0:N-1],q,'o','LineWidth',2,'MarkerSize',8);
xlabel('X (position)');ylabel('wave sum f(x)');

for n=2:length(q)
    figure(1);subplot(311);grid on;box on %real (cosine) part
    hold on;plot(X,cosine(X,real(fk(n)),n-1,N),lnsp{n},'LineWidth',2)
    title(['wave forms for k=',int2str(n)]);ylabel('real(f_k)');
     figure(1);subplot(312);grid on;box on %imaginary (sine) part
     hold on;plot(X,sine(X,imag(fk(n)),n-1,N),lnsp{n},'LineWidth',2)
     ylabel('imag(f_k)');
     wavesum=wavesum+cosine(X,real(fk(n)),n-1,N)-sine(X,imag(fk(n)),n-1,N);
     figure(1);subplot(313);hold off;grid on;box on
     plot([0:N-1],q,'o',X,wavesum+fk(1),'-k','LineWidth',2,'MarkerSize',8);
     xlabel('X (position)');ylabel('wave sum f(x)');
end