% clc
% clear all
% close all

%   强迫激励源


% Esource=exp( -((n-800)/100).^2);

% Ess=sin(2*pi*5e9*n*1/5e9/200);
% Es=10*exp( -((n-8000)/100).^2).*sin(2*pi*5e9*(n-800)*1e-9/3000);
% % subplot(3,2,1)
% 
% figure
% Es=exp( -((n-51)/17).^2);
% plot(n,Es)
% figure
% N=2^nextpow2(length(Vin));
% 
% Y=fft(Es,N);
% figure plot(abs(Y(1:N/2+1)/N))


N=2^nextpow2(length(Vin));
YYY=fft(Vin,N);
figure
plot(abs(YYY(1:N/2+1)))

