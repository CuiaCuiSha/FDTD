clc
clear all
close all


%   强迫激励源
% 

% % % 
% % % g = inline('t^2');
% % % argnames(g);
% % % a = 1;
% % % b = 2;
% % % c = 30;
% % % f = @(x,y) a*x.^2 + b*x + c*y;
% % % a=[1,1];
% % % b=[1,2];
% % % f(1, b)

n=1:1733;

% Esource=exp( -((n-800)/100).^2);

% Ess=sin(2*pi*5e9*n*1/5e9/200);
% Es=10*exp( -((n-8000)/100).^2).*sin(2*pi*5e9*(n-800)*1e-9/3000);
% % subplot(3,2,1)


Es=10*exp( -((n-300)/100).^2);
plot(n,Es)

% Fs=fft(Ess);
% Fsab=abs(Fs);
% figure

plot(linspace(0,5e9,length(Fsab)),Fsab)
% figure

% plot(exp(-(pi*1000*).^2))



