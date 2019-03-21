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
fmax=5e9;     %   3 mm

dx=1/20*c/fmax;     X_num=fix(X_long/dx);
dy=1/20*c/fmax;     Y_num=fix(Y_long/dy);
dz=1/20*c/fmax;     Z_num=fix(Z_long/dz);

%稳定性条件，dxdydz相等
dt=1/( c*sqrt(3/(dx)^2) );
TimeLong=fix(3* max([X_long,Y_long,Z_long]) /c/dt );
%----------------------------------%
%   dt 应当符合稳定性条件
%   dt leq 1/c*sqrt（1/dx2+1/dy2+1/dz2）
%   5 Ghz ，dt less 0.0021 s
%----------------------------------%
%   激励源设定，
t_max=1/2/fmax;
t_decay=fix(t_max/dt)
t0=4*t_decay;
t_source=6*t_decay;

%%   场的初始化
%----------------------------------%
Ex=zeros(X_num,Z_num+1);
Ez=zeros(X_num+1,Z_num);
Hy=zeros(X_num,Z_num);

Ext=Ex;Ezt=Ez;Hyt=Hy;

%%   计算
for t=1:3800
    %%   source
    if t<= t_source
        
        Ex(1:X_num,50)=10*exp( -((t-t0)/t_decay).^2);
        
    end
    
    %% 蛙跳算法
    % 保存上一时刻，和上上一时刻的场用作边界条件
    Extt=Ext;
    Extt=Ezt;
    Hytt=Hyt;
    Ext=Ex;
    Ezt=Ez;
    Hyt=Hy;
    
 
    
    
    Hy=Hy+dt/mu*(  ( Ez(2:X_num+1,:)-Ez(1:X_num,:) )./dx + ( Ex(:,1:Z_num)-Ex(:,2:Z_num+1) )./dz) ;
    Ex(:,2:Z_num)=Ex(:,2:Z_num)+...
        dt/epsilon*( ( Hy(:,1:Z_num-1)-Hy(:,2:Z_num) )/dz  );
    Ez(2:X_num,:)=Ez(2:X_num,:)+...
        dt/epsilon*( ( Hy(2:X_num,:)-Hy(1:X_num-1,:) )/dx );
    
    %% 边界条件
    %   上下壁板
    Ez(1,:)=0;
    Ez(X_num+1,:)=0;
    %   挡板
    Ex(1:fix(X_num/3),fix(Z_num/2))=0;
    Ex(fix(2*X_num/3):X_num,fix(Z_num/2))=0;
    %
    % 一阶Mur吸收条件 %吸收条件中真正起作用的其实是Ex，Ez的吸收
    Ex(:,1)=Ext(:,2)+( (c*dt-dz)/(c*dt+dz)*( Ex(:,2)-Ext(:,1)) );
    Ez(:,1)=Ezt(:,2)+( (c*dt-dz)/(c*dt+dz)*( Ez(:,2)-Ezt(:,1)) );
    %     Hy(:,1)=Hyt(:,2)+( (c*dt-dz)/(c*dt+dz)*( Hy(:,2)-Hyt(:,1)) );
    %
    
    Ex(:,Z_num+1)=Ext(:,Z_num)+( (c*dt-dz)/(c*dt+dz)*( Ex(:,Z_num)-Ext(:,Z_num+1)) );
    Ez(:,Z_num)=Ezt(:,Z_num-1)+( (c*dt-dz)/(c*dt+dz)*( Ez(:,Z_num-1)-Ezt(:,Z_num)) );
    %     Hy(:,Z_num)=Hyt(:,Z_num-1)+( (c*dt-dz)/(c*dt+dz)*( Hy(:,Z_num-1)-Hyt(:,Z_num)) );
     
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
    
    %% 绘动图
    PEx=Ex(1:X_num,1:Z_num);
    PEz=Ez(1:X_num,1:Z_num);
    Eabs=sqrt(PEx.^2+PEz.^2);
    [xx,yy]=meshgrid(1:Z_num,1:X_num);
    
    
    mesh(xx,yy,Eabs)
    view(0,90)
    pause(0.00000001)
    
    
end

%% 分离入射电压反射电压
%   方法比较笨。计算的是，Vin的max所在，作为入射电压中点，前后个取

[value,VinMid]=max(Verf_L);
[value,VinrefMid]=min(Verf_L);

%% 后处理
figure(2)
subplot(3,1,1)
plot(Verf_L);hold on;plot(Verf_R);title(' 左右电压波形 ');
legend('左侧挡板电压波形','右侧挡板电压波形')
subplot(3,1,2)
plot(Verf_L);title(' 左侧挡板电压波形 ')
subplot(3,1,3)
plot(Verf_R);title(' 右侧挡板电压波形 ')

figure(3)
subplot(3,1,1)
plot(Verf_L(VinMid-50:VinMid+50));

subplot(3,1,2)
plot(Verf_L(VinrefMid-50:VinrefMid+50));



