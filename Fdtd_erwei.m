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
epsilon=8.854e-12;


%----------------------------------%
%   X is long /m along axis x Y is wide /m along axis y
X_long=1;
Y_long=1;
Z_long=1;

Freq_max=5e9;     %   3 mm

dx=1/20*c/Freq_max;
dy=1/20*c/Freq_max;
% dz=1/20*c/Freq_max;

X_num=fix(X_long/dx);
Y_num=fix(Y_long/dy);
% Z_num=fix(Z_long/dz);



%稳定性条件，dxdydz相等
dt=1/( c*sqrt(3/(dx)^2) );
TimeLong=fix(3* max([X_long,Y_long,Z_long]) /c/dt );
%----------------------------------%
%   dt 应当符合稳定性条件
%   dt leq 1/c*sqrt（1/dx2+1/dy2+1/dz2）
%   5 Ghz ，dt less 0.0021 s
%----------------------------------%



%----------------------------------%
%   初始化
%----------------------------------%
Ex=zeros(X_num,Y_num+1,TimeLong+1);
Ey=zeros(X_num+1,Y_num,TimeLong+1);
Ez=zeros(X_num+1,Y_num+1,TimeLong+1);

Hx=zeros(X_num+1,Y_num,TimeLong+1);
Hy=zeros(X_num,Y_num+1,TimeLong+1);
Hz=zeros(X_num,Y_num,TimeLong+1);


%----------------------------------%
%   计算
%----------------------------------%

% % %   引入z_flag,三维--二维的
% % z_flag=1;
for t=1:TimeLong
    
    
    %蛙跳算法
    
    for jj=2:Y_num
        
        for ii=2:X_num
            
            Hx(ii,jj,kk,t+1)=Hx(ii,jj,kk,t)+...
                dt/mu*( ( Ey(ii,jj,kk+1,t)-Ey(ii,jj,kk,t) )/dz+...
                ( Ez(ii,jj,kk,t)-Ez(ii,jj+1,kk,t) )/dy  );
            
            Hy(ii,jj,kk,t+1)=Hy(ii,jj,kk,t)+...
                dt/mu*( ( Ez(ii+1,jj,kk,t)-Ez(ii,jj,kk,t) )/dx +...
                ( Ex(ii,jj,kk,t)-Ex(ii,jj,kk+1,t) )/dz  );
            
            
            Hz(ii,jj,kk,t+1)=Hz(ii,jj,kk,t)+...
                dt/mu*( ( Ex(ii,jj+1,kk,t)-Ex(ii,jj,kk,t)  )/dy+...
                ( Ey(ii,jj,kk,t)-Ey(ii+1,jj,kk,t)  )/dx  );
            
            
            Ex(ii,jj,kk,t+1)=Ex(ii,jj,kk,t)+...
                dt/epsilon*( ( Hz(ii,jj,kk,t+1)-Hz(ii,jj-1,kk,t+1) )/dy+...
                ( Hy(ii,jj,kk-1,t+1)-Hy(ii,jj,kk,t+1) )/dz  );
            
            Ey(ii,jj,kk,t+1)=Ey(ii,jj,kk,t)+...
                dt/epsilon*( ( Hx(ii,jj,kk,t+1)-Hx(ii,jj,kk-1,t+1) )/dz+...
                ( Hz(ii-1,jj,kk,t+1)-Hz(ii,jj,kk,t+1) )/dx  );
            
            Ez(ii,jj,kk,t+1)=Ez(ii,jj,kk,t)+...
                dt/epsilon*( ( Hy(ii,jj,kk,t+1)-Hy(ii-1,jj,kk,t+1) )/dx+...
                ( Hx(ii,jj-1,kk,t+1)-Hx(ii,jj,kk,t+1) )/dy  );
            
            
        end
        
    end
    
    
    %边界条件

   
    
    
    
    
    
end




