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

dx=1/20*c/Freq_max
dy=1/20*c/Freq_max
dz=1/20*c/Freq_max ;

X_num=fix(X_long/dx);
Y_num=fix(Y_long/dy);
Z_num=fix(Z_long/dz);
Z_num=1;

TimeLong=1000;     
%----------------------------------%
%   dt 应当符合稳定性条件
%   dt leq 1/sqrt（1/dx2+1/dy2+1/dz2）
%   5 Ghz ，dt less 0.0021 s 
%----------------------------------%
dt=0.0015;



%----------------------------------%
%   初始化
%----------------------------------%
Ex=zeros(X_num+1,Y_num+1,Z_num+1,TimeLong+1);
Ey=zeros(X_num+1,Y_num+1,Z_num+1,TimeLong+1);
Ez=zeros(X_num+1,Y_num+1,Z_num+1,TimeLong+1);

Hx=zeros(X_num,Y_num,Z_num,TimeLong+1);
Hy=zeros(X_num,Y_num,Z_num,TimeLong+1);
Hz=zeros(X_num,Y_num,Z_num,TimeLong+1);



%----------------------------------%
%   计算
%----------------------------------%

%   引入z_flag,三维--二维的
z_flag=1;
for t=1:TimeLong
    
    
    
    %蛙跳
    %为了保证四维数组有意义，t从1开始取值
    %   tips 磁场H，-1/2不计，+1/2，认为+1   用ceil向上取整
    %        电场E，+1/2不计                用floor向下取整
    
    for kk=1:Z_num
    
        for jj=1:Y_num
        
            for ii=1:X_num
                
                Hx(,,,,t+1)=Hx(,,,t)+...
                    dt/mu*( ( Ey(,,,t)-Ey(,,,t) )/dz+...
                            ( Ez(,,,t)-Ez(,,,t) )/dy  );
                
                Hy(,,,t+1)=Hy(,,,t)+...
                    dt/mu*( ( Ez(,,,t)-Ez(,,,t) )/dx +...
                            ( Ex(floor(ii+0.5),jj,kk,t)-Ex(floor(ii+0.5),jj,kk+1,t) )/dz  );
                
                
                Hz(ceil(ii+0.5),ceil(jj+0.5),kk,t+1)=Hz(ceil(ii+0.5),ceil(jj+0.5),kk,t)+...
                    dt/mu*( ( Ex(floor(ii+0.5),jj+1,kk,t)-Ex(floor(ii+0.5),jj,kk,t)  )/dy+...
                            ( Ey(ii,floor(jj+0.5),kk,t)-Ey(ii-1,floor(jj+0.5),kk,t)  )/dx  );          

                     
                Ex(floor(ii+0.5),jj,kk,t+1)=Ex(floor(ii+0.5),jj,kk,t)+...
                    dt/epsilon*( ( Hz(,,kk,t+1)-Hy(,,kk,t+1) )/dy+...
                                 ( Hz(,jj,,t+1)-Hz(,jj,,t+1) )/dz  );
                
                Ey(,,,t+1)=Ey(,,,t)+...
                    dt/epsilon*( ( Hx(,,,t+1)-Hx(,,,t+1) )/dz+...
                                 ( Hz(,,,t+1)-Hz(,,,t+1) )/dx  );
                
                Ez(,,,t+1)=Ez(,,,t)+...
                    dt/epsilon*( ( Hy(,,,t+1)-Hy(,,,t+1) )/dx+...
                                 ( Hx(,,,t+1)-Hx(,,,t+1) )/dy  );
                
% % %                         
% % %                 Hy(ceil(ii+0.5),jj,ceil(kk+0.5),t+1)=Hy(ceil(ii+0.5),jj,ceil(kk+0.5),t)+...
% % %                     dt/mu*( ( Ez(ii+1,jj,floor(kk+0.5),t)-Ez(ii,jj,floor(kk+0.5),t) )/dx +...
% % %                             ( Ex(floor(ii+0.5),jj,kk,t)-Ex(floor(ii+0.5),jj,kk+1,t) )/dz  );
% % %                 
% % %                 
% % %                 Hz(ceil(ii+0.5),ceil(jj+0.5),kk,t+1)=Hz(ceil(ii+0.5),ceil(jj+0.5),kk,t)+...
% % %                     dt/mu*( ( Ex(floor(ii+0.5),jj+1,kk,t)-Ex(floor(ii+0.5),jj,kk,t)  )/dy+...
% % %                             ( Ey(ii,floor(jj+0.5),kk,t)-Ey(ii-1,floor(jj+0.5),kk,t)  )/dx  );                
                
            end
            
        end
        
    end
    
    
    
    
    %边界条件
    
    
    
end




