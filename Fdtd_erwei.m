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
                dt/epsilon*( ( Hz(ii,jj,kk,t+1)-Hy(ii,jj-1,kk,t+1) )/dy+...
                ( Hz(ii,jj,kk-1,t+1)-Hz(ii,jj,kk,t+1) )/dz  );
            
            Ey(ii,jj,kk,t+1)=Ey(ii,jj,kk,t)+...
                dt/epsilon*( ( Hx(ii,jj,kk,t+1)-Hx(ii,jj,kk-1,t+1) )/dz+...
                ( Hz(ii-1,jj,kk,t+1)-Hz(ii,jj,kk,t+1) )/dx  );
            
            Ez(ii,jj,kk,t+1)=Ez(ii,jj,kk,t)+...
                dt/epsilon*( ( Hy(ii,jj,kk,t+1)-Hy(ii-1,jj,kk,t+1) )/dx+...
                ( Hx(ii,jj-1,kk,t+1)-Hx(ii,jj,kk,t+1) )/dy  );
            
            
        end
        
    end
    
    
    %边界条件
    %三阶廖氏吸收边界。
    if t>2
        Ex(1,:,t+1)=8/9*Ex(1,:,t)+4/9*Ex(2,:,t)-3/8*Ex(3,:,t)-3*Ex(2,:,t-1)-...
            1/8*Ex(1,:,t-2)+3/4*Ex(2,:,t-2)+3/8*Ex(3,:,t-2);
        
        Ey(1,:,t+1)=8/9*Ey(1,:,t)+4/9*Ey(2,:,t)-3/8*Ey(3,:,t)-3*Ey(2,:,t-1)-...
            1/8*Ey(1,t-2)+3/4*Ey(2,:,t-2)+3/8*Ey(3,:,t-2);
        
        Ez(1,:,t+1)=8/9*Ez(1,:,t)+4/9*Ez(2,:,t)-3/8*Ez(3,:,t)-3*Ez(2,:,t-1)-...
            1/8*Ez(1,:,t-2)+3/4*Ez(2,:,t-2)+3/8*Ez(3,:,t-2);
        
        Hx(1,:,t+1)=8/9*Hx(1,:,t)+4/9*Hx(2,:,t)-3/8*Hx(3,:,t)-3*Hx(2,:,t-1)-...
            1/8*Hx(1,:,t-2)+3/4*Hx(2,:,t-2)+3/8*Hx(3,:,t-2);
        
        Hy(1,:,t+1)=8/9*Hy(1,:,t)+4/9*Hy(2,:,t)-3/8*Hy(3,:,t)-3*Hy(2,:,t-1)-...
            1/8*Hy(1,:,t-2)+3/4*Hy(2,:,t-2)+3/8*Hy(3,:,t-2);
        
        Hz(1,:,t+1)=8/9*Hz(1,:,t)+4/9*Hz(2,:,t)-3/8*Hz(3,:,t)-3*Hz(2,:,t-1)-...
            1/8*Hz(1,:,t-2)+3/4*Hz(2,:,t-2)+3/8*Hz(3,:,t-2);
    elseif t==2
        
        Ex(1,:,t+1)=8/9*Ex(1,:,t)+4/9*Ex(2,:,t)-3/8*Ex(3,:,t)-3*Ex(2,:,t-1);
        Ey(1,:,t+1)=8/9*Ey(1,:,t)+4/9*Ey(2,:,t)-3/8*Ey(3,:,t)-3*Ey(2,:,t-1);
        Ez(1,:,t+1)=8/9*Ez(1,:,t)+4/9*Ez(2,:,t)-3/8*Ez(3,:,t)-3*Ez(2,:,t-1);
        Hx(1,:,t+1)=8/9*Hx(1,:,t)+4/9*Hx(2,:,t)-3/8*Hx(3,:,t)-3*Hx(2,:,t-1);
        Hy(1,:,t+1)=8/9*Hy(1,:,t)+4/9*Hy(2,:,t)-3/8*Hy(3,:,t)-3*Hy(2,:,t-1);
        Hz(1,:,t+1)=8/9*Hz(1,:,t)+4/9*Hz(2,:,t)-3/8*Hz(3,:,t)-3*Hz(2,:,t-1);
    else
        
        Ex(1,:,t+1)=8/9*Ex(1,:,t)+4/9*Ex(2,:,t)-3/8*Ex(3,:,t);
        Ey(1,:,t+1)=8/9*Ey(1,:,t)+4/9*Ey(2,:,t)-3/8*Ey(3,:,t);
        Ez(1,:,t+1)=8/9*Ez(1,:,t)+4/9*Ez(2,:,t)-3/8*Ez(3,:,t);
        Hx(1,:,t+1)=8/9*Hx(1,:,t)+4/9*Hx(2,:,t)-3/8*Hx(3,:,t);
        Hy(1,:,t+1)=8/9*Hy(1,:,t)+4/9*Hy(2,:,t)-3/8*Hy(3,:,t);
        Hz(1,:,t+1)=8/9*Hz(1,:,t)+4/9*Hz(2,:,t)-3/8*Hz(3,:,t);
        
    end
        
    
    
end




