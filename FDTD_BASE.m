clc
clear all
close all


%----------------------------------%
%   3-7 Ghz     5 Ghz
%   
%   data:2018.10.17
%----------------------------------%


%----------------------------------%
%   常数项设定
%----------------------------------%

c=2.998e8;
mu=4*pi*1e-7;
vareps=8.854e-12;


%----------------------------------%
%   X is long /m along axis x
%   Y is wide /m along axis y
X_long=1;
Y_long=1;
Z_long=1;

freq_max=5e9;     %   3 mm

dx=1/20*c/freq_max 
dy=1/20*c/freq_max 
dz=1/20*c/freq_max ;

X_num=fix(X_long/dx);
Y_num=fix(Y_long/dy);
Z_num=fix(Z_long/dz);
Z_num=1;
%----------------------------------%
%   
%----------------------------------%


%----------------------------------%
%   初始化
%----------------------------------%
Ex=zeros(X_num,Y_num,Z_num);
Ey=zeros(X_num,Y_num,Z_num);
Ez=zeros(X_num,Y_num,Z_num);

Hx=zeros(X_num,Y_num,Z_num);
Hy=zeros(X_num,Y_num,Z_num);
Hz=zeros(X_num,Y_num,Z_num);

TimeLong=1000;


for t=0:TimeLong







end

