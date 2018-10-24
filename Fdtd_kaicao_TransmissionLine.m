%%   clear
%----------------------------------%
%   3-7 Ghz     5 Ghz
%
%   data:2018.10.17
%----------------------------------%
clc;clear all;close all
%----------------------------------%
%%   常数项，固定量设定 
%----------------------------------%

c=2.998e8;
mu=4*pi*1e-7;
epsilon=8.854e-12;

%----------------------------------%
%   X is long /m along axis x Y is wide /m along axis y
X_long=0.3;
Y_long=1;
Z_long=1;

Freq_max=5e9;     %   3 mm

dx=1/20*c/Freq_max;     X_num=fix(X_long/dx);
dy=1/20*c/Freq_max;     Y_num=fix(Y_long/dy);
dz=1/20*c/Freq_max;     Z_num=fix(Z_long/dz);

%稳定性条件，dxdydz相等
dt=1/( c*sqrt(3/(dx)^2) );
TimeLong=fix(3* max([X_long,Y_long,Z_long]) /c/dt );
%----------------------------------%
%   dt 应当符合稳定性条件
%   dt leq 1/c*sqrt（1/dx2+1/dy2+1/dz2）
%   5 Ghz ，dt less 0.0021 s
%----------------------------------%

%----------------------------------%
%%   场的初始化
%----------------------------------%
Ex=zeros(X_num,Z_num+1);
% Ey=zeros(X_num+1,Z_num+1,TimeLong+1);
Ez=zeros(X_num+1,Z_num);
% Hx=zeros(X_num+1,Z_num,TimeLong+1);
Hy=zeros(X_num,Z_num);
% Hz=zeros(X_num,Z_num+1,TimeLong+1);

%%   计算
for t=1:600
    %%   source
    if t<= 60
        
        
%         Ex(1:X_num,10:40)=0;
%         Ez(1:X_num,10:40)=0;
%         Hy(1:X_num,10:40)=0;
        Ex(1:X_num,50)=10*exp( -((t-30)/10).^2);
        
    end
%     if t==200
%         Ex(1:X_num,10:100)=0;
%         Ez(1:X_num,10:100)=0;
%         Hy(1:X_num,10:100)=0;
%     end
    
    %% 蛙跳算法
    % 保存上一时刻的场用作边界条件
    Ext=Ex;
    Ezt=Ez;
    Hyt=Hy;
    
    for kk=2:Z_num
        for ii=2:X_num
            Hy(ii,kk)=Hy(ii,kk)+...
                dt/mu*( ( Ez(ii+1,kk)-Ez(ii,kk) )/dx +...
                ( Ex(ii,kk)-Ex(ii,kk+1) )/dz  );
        end
    end
    for kk=2:Z_num
        Ex(:,kk)=Ex(:,kk)+...
            dt/epsilon*( ( Hy(:,kk-1)-Hy(:,kk) )/dz  );
    end
    
    for ii=2:X_num
        Ez(ii,:)=Ez(ii,:)+...
            dt/epsilon*( ( Hy(ii,:)-Hy(ii-1,:) )/dx );
    end
    
    %             Hx(ii,jj,kk,t+1)=Hx(ii,jj,kk,t)+...
    %                 dt/mu*( ( Ey(ii,jj,kk+1,t)-Ey(ii,jj,kk,t) )/dz+...
    %                 ( Ez(ii,jj,kk,t)-Ez(ii,jj+1,kk,t) )/dy  );
    %             Hz(ii,jj,kk,t+1)=Hz(ii,jj,kk,t)+...
    %                 dt/mu*( ( Ex(ii,jj+1,kk,t)-Ex(ii,jj,kk,t)  )/dy+...
    %                 ( Ey(ii,jj,kk,t)-Ey(ii+1,jj,kk,t)  )/dx  );
    %             Ey(ii,jj,kk,t+1)=Ey(ii,jj,kk,t)+...
    %                 dt/epsilon*( ( Hx(ii,jj,kk,t+1)-Hx(ii,jj,kk-1,t+1) )/dz+...
    %                 ( Hz(ii-1,jj,kk,t+1)-Hz(ii,jj,kk,t+1) )/dx  );
    
    %% 边界条件
%     Ez(1,:)=0;
%     Ez(X_num+1,:)=0;
%     Ex(:,1)=0;
%     Ex(:,Z_num+1)=0;  
    %   上下壁板
    Ez(1:3,:)=0;
    Ez(X_num-2:X_num,:)=0;
    
    %   挡板
    Ex(1:fix(X_num/3),fix(Z_num/2))=0;
    Ex(fix(2*X_num/3):X_num,fix(Z_num/2))=0;
%     Ez(1:fix(X_num/3),fix(Z_num/2):fix(Z_num/2)+3)=0;
%     Ez(fix(2*X_num/3):X_num,fix(Z_num/2):fix(Z_num/2)+3)=0;
    
    %% 一阶Mur吸收条件
    Ex(:,1)=Ext(:,2)+( (c*dt-dx)/(c*dt+dx)*(Ex(:,2)-Ext(:,1)) );
    Ez(:,1)=Ezt(:,2)+( (c*dt-dz)/(c*dt+dz)*(Ez(:,2)-Ezt(:,1)) );
    Hy(:,1)=Hyt(:,2)+( (c*dt-dy)/(c*dt+dy)*(Hy(:,2)-Hyt(:,1)) );
    
    %% 二阶Mur吸收条件
%     Ex(:,1)=-Ext(:,2)+(c*dt-dx)/(c*dt+dx)*(Ex(:,2)-Ext(:,1)) +2*dx/(c*dt+dx)*(Ext(:,1)+Ext(:,2))
    
    %% 廖氏边界条件 已屏蔽
    %
    %     if t>-1
    %         %                 Ex(1,:,t+1)=8/9*Ex(1,:,t)+4/9*Ex(2,:,t)-3/8*Ex(3,:,t)-3*Ex(2,:,t-1)-...
    %         %                     1/8*Ex(1,:,t-2)+3/4*Ex(2,:,t-2)+3/8*Ex(3,:,t-2);
    %         Ex(:,1,t+1)=8/9*Ex(:,1,t)+4/9*Ex(:,2,t)-3/8*Ex(:,3,t)-3*Ex(:,2,t-1)-...
    %             1/8*Ex(:,1,t-2)+3/4*Ex(:,2,t-2)+3/8*Ex(:,3,t-2);
    %         %
    %         %         Ey(1,:,t+1)=8/9*Ey(1,:,t)+4/9*Ey(2,:,t)-3/8*Ey(3,:,t)-3*Ey(2,:,t-1)-...
    %         %             1/8*Ey(1,t-2)+3/4*Ey(2,:,t-2)+3/8*Ey(3,:,t-2);
    %
    %         %                 Ez(1,:,t+1)=8/9*Ez(1,:,t)+4/9*Ez(2,:,t)-3/8*Ez(3,:,t)-3*Ez(2,:,t-1)-...
    %         %                     1/8*Ez(1,:,t-2)+3/4*Ez(2,:,t-2)+3/8*Ez(3,:,t-2);
    %         Ez(:,1,t+1) =8/9*Ez(:,1,t)+4/9*Ez(:,2,t)-3/8*Ez(:,3,t)-3*Ez(:,2,t-1)-...
    %             1/8*Ez(:,1,t-2)+3/4*Ez(:,2,t-2)+3/8*Ez(:,3,t-2);
    %         %         Hx(1,:,t+1)=8/9*Hx(1,:,t)+4/9*Hx(2,:,t)-3/8*Hx(3,:,t)-3*Hx(2,:,t-1)-...
    %         %             1/8*Hx(1,:,t-2)+3/4*Hx(2,:,t-2)+3/8*Hx(3,:,t-2);
    %
    %         %                 Hy(1,:,t+1)=8/9*Hy(1,:,t)+4/9*Hy(2,:,t)-3/8*Hy(3,:,t)-3*Hy(2,:,t-1)-...
    %         %                     1/8*Hy(1,:,t-2)+3/4*Hy(2,:,t-2)+3/8*Hy(3,:,t-2);
    %
    %         Hy(:,1,t+1)=8/9*Hy(:,1,t)+4/9*Hy(:,2,t)-3/8*Hy(:,3,t)-3*Hy(:,2,t-1)-...
    %             1/8*Hy(:,1,t-2)+3/4*Hy(:,2,t-2)+3/8*Hy(:,3,t-2);
    %         %
    %         %         Hz(1,:,t+1)=8/9*Hz(1,:,t)+4/9*Hz(2,:,t)-3/8*Hz(3,:,t)-3*Hz(2,:,t-1)-...
    %         %             1/8*Hz(1,:,t-2)+3/4*Hz(2,:,t-2)+3/8*Hz(3,:,t-2);
    %
    %
    %
    %         %     elseif t==2
    % %
    % %         %                 Ex(1,:,t+1)=8/9*Ex(1,:,t)+4/9*Ex(2,:,t)-3/8*Ex(3,:,t)-3*Ex(2,:,t-1);
    % %         Ex(:,1,t+1)=8/9*Ex(:,1,t)+4/9*Ex(:,2,t)-3/8*Ex(:,3,t)-3*Ex(:,2,t-1);
    % %         %         Ey(1,:,t+1)=8/9*Ey(1,:,t)+4/9*Ey(2,:,t)-3/8*Ey(3,:,t)-3*Ey(2,:,t-1);
    % %         %                 Ez(1,:,t+1)=8/9*Ez(1,:,t)+4/9*Ez(2,:,t)-3/8*Ez(3,:,t)-3*Ez(2,:,t-1);
    % %         Ez(:,1,t+1)=8/9*Ez(:,1,t)+4/9*Ez(:,2,t)-3/8*Ez(:,3,t)-3*Ez(:,2,t-1);
    % %         %         Hx(1,:,t+1)=8/9*Hx(1,:,t)+4/9*Hx(2,:,t)-3/8*Hx(3,:,t)-3*Hx(2,:,t-1);
    % %         %                 Hy(1,:,t+1)=8/9*Hy(1,:,t)+4/9*Hy(2,:,t)-3/8*Hy(3,:,t)-3*Hy(2,:,t-1);
    % %         Hy(:,1,t+1)=8/9*Hy(:,1,t)+4/9*Hy(:,2,t)-3/8*Hy(:,3,t)-3*Hy(:,2,t-1);
    % %         %         Hz(1,:,t+1)=8/9*Hz(1,:,t)+4/9*Hz(2,:,t)-3/8*Hz(3,:,t)-3*Hz(2,:,t-1);
    % %     else
    % %
    % %         %                 Ex(1,:,t+1)=8/9*Ex(1,:,t)+4/9*Ex(2,:,t)-3/8*Ex(3,:,t);
    % %         Ex(:,1,t+1)=8/9*Ex(:,1,t)+4/9*Ex(:,2,t)-3/8*Ex(:,3,t);
    % %         %         Ey(1,:,t+1)=8/9*Ey(1,:,t)+4/9*Ey(2,:,t)-3/8*Ey(3,:,t);
    % %         %                 Ez(1,:,t+1)=8/9*Ez(1,:,t)+4/9*Ez(2,:,t)-3/8*Ez(3,:,t);
    % %         Ez(:,1,t+1)=8/9*Ez(:,1,t)+4/9*Ez(:,2,t)-3/8*Ez(:,3,t);
    % %         %         Hx(1,:,t+1)=8/9*Hx(1,:,t)+4/9*Hx(2,:,t)-3/8*Hx(3,:,t);
    % %         %                 Hy(1,:,t+1)=8/9*Hy(1,:,t)+4/9*Hy(2,:,t)-3/8*Hy(3,:,t);
    % %         Hy(:,1,t+1)=8/9*Hy(:,1,t)+4/9*Hy(:,2,t)-3/8*Hy(:,3,t);
    % %         %         Hz(1,:,t+1)=8/9*Hz(1,:,t)+4/9*Hz(2,:,t)-3/8*Hz(3,:,t);
    % %
    %     end
    
    %% 记录电压波形
    Verf_L(t)=sum( Ext(4:X_num-3,fix(Z_num/4)) )*dx;
    Verf_R(t)=sum( Ext(4:X_num-3,fix(Z_num*3/4)) )*dx;
    Verf_BAN(t)=sum( Ext(4:X_num-3,fix(Z_num/2)) )*dx;
    
    %% 绘图
    PEx=Ex(1:X_num,1:Z_num);
    PEz=Ez(1:X_num,1:Z_num);
    Eabs=sqrt(PEx.^2+PEz.^2);
    [xx,yy]=meshgrid(1:Z_num,1:X_num);
    
    
    mesh(xx,yy,PEx.^2)
    view(0,90)
    pause(0.00000001)
    
end

figure 
plot(Verf_L);
hold on
plot(Verf_R);
hold on
plot(Verf_BAN);
legend('Left','Righr','ban')



