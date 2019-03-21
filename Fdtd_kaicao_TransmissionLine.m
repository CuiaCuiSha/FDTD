
%----------------------------------%
%   3-7 Ghz     5 Ghz
%
%   data:2018.10.17
%----------------------------------%

clc;clear all;close all
feature('DefaultCharacterSet','UTF-8');
%----------------------------------%
%%  常数项，固定量设定
%----------------------------------%
c=2.998e8;
mu=4*pi*1e-7;
epsilon=8.854e-12;
%----------------------------------%
%   X is long /m along axis x Y is wide /m along axis y
X_long=1;
Y_long=1;
Z_long=1;
fmax=3e9;     %   3 mm

dx=1/20*c/fmax;     X_num=fix(X_long/dx);
dy=1/20*c/fmax;     Y_num=fix(Y_long/dy);
dz=1/20*c/fmax;     Z_num=fix(Z_long/dz);

%稳定性条件，dxdydz相等
dt=1/( 2*c*sqrt(3/(dx)^2) );
TimeLong=fix(3* max([X_long,Y_long,Z_long]) /c/dt );
%----------------------------------%
%   dt 应当符合稳定性条件
%   dt leq 1/c*sqrt（1/dx2+1/dy2+1/dz2）
%   5 Ghz ，dt less 0.0021 s
%----------------------------------%
%   激励源设定，
t_max=1/2/fmax;
t_decay=fix(t_max/dt);
t0=4*t_decay;
t_source=6*t_decay;

%%  第一次-场的初始化
%----------------------------------%
Ex=zeros(X_num,Z_num+1);
Ez=zeros(X_num+1,Z_num);
Hy=zeros(X_num,Z_num);
% TRME1=zeros()



%预分配提高速度用
Vref_Ltotal=zeros(1,TimeLong);
Vref_Rtotal=zeros(1,TimeLong);
Vref_BAN=zeros(1,TimeLong);
Ext=Ex;Ezt=Ez;Hyt=Hy;

%%  第一次-计算
for t=1:TimeLong
    %%   source
    if t< t_source
        Hy(100,100)=1*exp( -((t-t0)/t_decay).^2);
    end
    if  t== 3*t_source
       close all
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
    %   左右壁板
    Ex(:,3)=0;
    Ex(:,Z_num-3)=0;
    %   有厚度挡板
%     Ex(1:fix(X_num/3),fix(Z_num/2)-2:fix(Z_num/2)+2)=0;
%     Ex(fix(2*X_num/3):X_num,fix(Z_num/2)-2:fix(Z_num/2)+2)=0;
    %   无厚度挡板
%     Ex(1:fix(X_num/3),fix(Z_num/2):fix(Z_num/2)+1)=0;
%     Ex(fix(2*X_num/3):X_num,fix(Z_num/2):fix(Z_num/2)+1)=0;
    %
    
    
    
    % 一阶Mur吸收条件 %吸收条件中真正起作用的其实是Ex，Ez的吸收
    Ex(:,1)=Ext(:,2)+( (c*dt-dz)/(c*dt+dz)*( Ex(:,2)-Ext(:,1)) );
    Ez(:,1)=Ezt(:,2)+( (c*dt-dz)/(c*dt+dz)*( Ez(:,2)-Ezt(:,1)) );
    %     Hy(:,1)=Hyt(:,2)+( (c*dt-dz)/(c*dt+dz)*( Hy(:,2)-Hyt(:,1)) );
    Ex(:,Z_num+1)=Ext(:,Z_num)+( (c*dt-dz)/(c*dt+dz)*( Ex(:,Z_num)-Ext(:,Z_num+1)) );
    Ez(:,Z_num)=Ezt(:,Z_num-1)+( (c*dt-dz)/(c*dt+dz)*( Ez(:,Z_num-1)-Ezt(:,Z_num)) );
    %     Hy(:,Z_num)=Hyt(:,Z_num-1)+( (c*dt-dz)/(c*dt+dz)*( Hy(:,Z_num-1)-Hyt(:,Z_num)) );
    
    
    %% 记录电压波形

    Vref_Ltotal(t)=sum( Ext(4:X_num-3,fix(Z_num/4)) )*dx;
    Vref_Rtotal(t)=sum( Ext(4:X_num-3,fix(Z_num*3/4)) )*dx;
    Vref_BAN(t)=sum( Ext(4:X_num-3,fix(Z_num/2)) )*dx;
    %   TRM的设置
    TrmE_L1(:,:,t)=Ext(4:X_num-3,fix(Z_num/4));
    TrmE_R1(:,:,t)=Ext(4:X_num-3,fix(Z_num*3/4));
    %% 绘动图
    PEx=Ex(1:X_num,1:Z_num);
    PEz=Ez(1:X_num,1:Z_num);
     Eabs=sqrt(PEx.^2+PEz.^2);
    [xx,yy]=meshgrid(1:Z_num,1:X_num);
    mesh(xx,yy,Eabs)
%     view(0,90)
    pause(0.00000001)
%     
end

%%  第一次-后处理绘图
figure
subplot(3,1,1)
plot(Vref_Ltotal);hold on;plot(Vref_Rtotal);title(' 左右电压波形 ');
legend('左侧挡板电压波形','右侧挡板电压波形')
subplot(3,1,2)
plot(Vref_Ltotal);title(' 左侧挡板电压波形 ')
subplot(3,1,3)
plot(Vref_Rtotal);title(' 右侧挡板电压波形 ')
suptitle('有挡板')
% % % 
% % % %%  第二次计算-场的初始化
% % % %----------------------------------%
% % % Ex=zeros(X_num,Z_num+1);
% % % Ez=zeros(X_num+1,Z_num);
% % % Hy=zeros(X_num,Z_num);
% % % 
% % % Ext=Ex;Ezt=Ez;Hyt=Hy;
% % % 
% % % %%  第二次计算-计算
% % % for t=1:TimeLong
% % %     %%   source
% % %     if t<= t_source
% % %         
% % %         Ex(1:X_num,50)=10*exp( -((t-t0)/t_decay).^2);
% % %         
% % %     end
% % %     
% % %     %% 蛙跳算法
% % %     % 保存上一时刻，和上上一时刻的场用作边界条件
% % %     Extt=Ext;
% % %     Extt=Ezt;
% % %     Hytt=Hyt;
% % %     Ext=Ex;
% % %     Ezt=Ez;
% % %     Hyt=Hy;
% % %     
% % %     Hy=Hy+dt/mu*(  ( Ez(2:X_num+1,:)-Ez(1:X_num,:) )./dx + ( Ex(:,1:Z_num)-Ex(:,2:Z_num+1) )./dz) ;
% % %     Ex(:,2:Z_num)=Ex(:,2:Z_num)+...
% % %         dt/epsilon*( ( Hy(:,1:Z_num-1)-Hy(:,2:Z_num) )/dz  );
% % %     Ez(2:X_num,:)=Ez(2:X_num,:)+...
% % %         dt/epsilon*( ( Hy(2:X_num,:)-Hy(1:X_num-1,:) )/dx );
% % %     
% % %     %% 边界条件
% % %     %   上下壁板
% % %     Ez(1,:)=0;
% % %     Ez(X_num+1,:)=0;
% % %     %     %   挡板
% % %     %     Ex(1:fix(X_num/3),fix(Z_num/2))=0;
% % %     %     Ex(fix(2*X_num/3):X_num,fix(Z_num/2))=0;
% % %     
% % %     % 一阶Mur吸收条件 %吸收条件中真正起作用的其实是Ex，Ez的吸收
% % %     Ex(:,1)=Ext(:,2)+( (c*dt-dz)/(c*dt+dz)*( Ex(:,2)-Ext(:,1)) );
% % %     Ez(:,1)=Ezt(:,2)+( (c*dt-dz)/(c*dt+dz)*( Ez(:,2)-Ezt(:,1)) );
% % %     %     Hy(:,1)=Hyt(:,2)+( (c*dt-dz)/(c*dt+dz)*( Hy(:,2)-Hyt(:,1)) );
% % %     
% % %     Ex(:,Z_num+1)=Ext(:,Z_num)+( (c*dt-dz)/(c*dt+dz)*( Ex(:,Z_num)-Ext(:,Z_num+1)) );
% % %     Ez(:,Z_num)=Ezt(:,Z_num-1)+( (c*dt-dz)/(c*dt+dz)*( Ez(:,Z_num-1)-Ezt(:,Z_num)) );
% % %     %     Hy(:,Z_num)=Hyt(:,Z_num-1)+( (c*dt-dz)/(c*dt+dz)*( Hy(:,Z_num-1)-Hyt(:,Z_num)) );
% % %     
% % %     %% 记录电压波形
% % %     %   TRM的设置
% % %     TrmE_L2(:,:,t)=Ext(4:X_num-3,fix(Z_num/4));
% % %     TrmE_R2(:,:,t)=Ext(4:X_num-3,fix(Z_num*3/4));
% % %     %   L2 为不加挡板的情况，因此L2为入射电压
% % %     Vref_L2(t)=sum( Ext(4:X_num-3,fix(Z_num/4)) )*dx;
% % %     Vref_R2(t)=sum( Ext(4:X_num-3,fix(Z_num*3/4)) )*dx;
% % %     Vref_BAN2(t)=sum( Ext(4:X_num-3,fix(Z_num/2)) )*dx;
% % %     
% % %     %% 绘动图-已经全部屏蔽
% % %     %     PEx=Ex(1:X_num,1:Z_num);
% % %     %     PEz=Ez(1:X_num,1:Z_num);
% % %     %     Eabs=sqrt(PEx.^2+PEz.^2);
% % %     %     [xx,yy]=meshgrid(1:Z_num,1:X_num);
% % %     %
% % %     %
% % %     %     mesh(xx,yy,Eabs)
% % %     %     view(0,90)
% % %     %     pause(0.00000001)
% % %     
% % % end
% % % 
% % % %%  第二次后处理-没有参考意义-屏蔽
% % % % % figure
% % % % % subplot(3,1,1)
% % % % % plot(Vref_L2);hold on;plot(Vref_R2);title(' 左右电压波形2 ');
% % % % % legend('左侧挡板电压波形2','右侧挡板电压波形2')
% % % % % subplot(3,1,2)
% % % % % plot(Vref_L2);title(' 左侧挡板 入射电压 ')
% % % % % subplot(3,1,3)
% % % % % plot(Vref_R2);title(' 右侧挡板 投射电压 ')
% % % 
% % % %%  两次对比绘图 有挡板-无挡板
% % % Vref_Lre=Vref_Ltotal-Vref_L2;
% % % Vref_Lre=Vref_Ltotal-Vref_L2;
% % % figure
% % % subplot(3,1,1)
% % % plot(Vref_L2);hold on;title(' 左侧挡板 入射电压 ');
% % % subplot(3,1,2)
% % % plot(Vref_Lre);hold on;title(' 左侧挡板 反射电压 ');
% % % subplot(3,1,3)
% % % plot(Vref_Rtotal);hold on;title(' 右侧挡板 透射电压 ');
% % % suptitle('有挡板-无挡板')
% % % 
% % % %%  离散傅里叶变换-自己写-已经屏蔽
% % % % Vin=Vref_L2;
% % % % Vre=Vref_Lre;
% % % % Vtr=Vref_Rtotal;
% % % % aa=1733;
% % % % df=fmax/aa;
% % % %
% % % % for k=1:aa
% % % %     Vinf=0;
% % % %     Vref=0;Vtrf=0;
% % % %     for n=1:TimeLong
% % % %         Vinf=Vinf+Vin(n)*exp(-1i*2*pi*k*df*n*dt);
% % % %         Vref=Vref+Vre(n)*exp(-1i*2*pi*k*df*n*dt);
% % % %         Vtrf=Vtrf+Vtr(n)*exp(-1i*2*pi*k*df*n*dt);
% % % %     end
% % % %     Vinfft(k)=dt*Vinf;
% % % %     Vrefft(k)=dt*Vref;
% % % %     Vtrfft(k)=dt*Vtrf;
% % % % end
% % % %
% % % % figure
% % % % plot(0:df:fmax-df,abs(Vtrfft(:))./abs(Vinfft(:)))
% % % % hold on
% % % % plot(0:df:fmax-df,abs(Vrefft(:))./abs(Vinfft(:)))
% % % % hold on
% % % % plot(0:df:fmax-df,   (abs(Vrefft(:))./abs(Vinfft(:))).^2+ (abs(Vtrfft(:))./abs(Vinfft(:))).^2   )
% % % 
% % % %%  直接采用自带函数fft
% % % clear N
% % % Vin=Vref_L2;
% % % Vre=Vref_Lre;
% % % Vtr=Vref_Rtotal;
% % % 
% % % N=2^nextpow2(length(Vin));
% % % Vinfft=fft(Vin,N);
% % % Vrefft=fft(Vre,N);
% % % Vtrfft=fft(Vtr,N);
% % % figure
% % % plot(abs(Vinfft(1:N/2+1)));hold on
% % % plot(abs(Vrefft(1:N/2+1)));hold on
% % % plot(abs(Vtrfft(1:N/2+1)));hold on
% % % legend('Vinfft','Vrefft','Vtrfft')
% % % % plot( sqrt(abs(Vtrfft(1:N/2+1)).^2+abs(Vrefft(1:N/2+1)).^2));
% % % figure
% % % plot(abs(Vrefft(1:N/2+1))./abs(Vinfft(1:N/2+1)));hold on
% % % plot(abs(Vtrfft(1:N/2+1))./abs(Vinfft(1:N/2+1)));hold on
% % % legend('S11','S21')
% % % 
% % % %%  第三次TR过程的计算-场的初始化
% % % % close all
% % % % 做差
% % % Trm_L=TrmE_L2(:,:,:)-TrmE_L1(:,:,:);
% % % Trm_R=TrmE_R2(:,:,:)-TrmE_R1(:,:,:);
% % % 
% % % %%  第三次计算-场的初始化
% % % %----------------------------------%
% % % Ex=zeros(X_num,Z_num+1);
% % % Ez=zeros(X_num+1,Z_num);
% % % Hy=zeros(X_num,Z_num);
% % % % TRME1=zeros()
% % % Ext=Ex;Ezt=Ez;Hyt=Hy;
% % % figure
% % % 
% % % %%  第三次计算-
% % % for t=1:TimeLong
% % %     %%   source
% % %   
% % %         
% % %     Ex(4:X_num-3,fix(Z_num/4))=TrmE_L1(:,:,TimeLong+1-t);
% % %     Ex(4:X_num-3,fix(Z_num*3/4))=TrmE_R1(:,:,TimeLong+1-t);
% % %         
% % %    
% % % 
% % %     
% % %     %% 蛙跳算法
% % %     % 保存上一时刻，和上上一时刻的场用作边界条件
% % %     Extt=Ext;
% % %     Extt=Ezt;
% % %     Hytt=Hyt;
% % %     Ext=Ex;
% % %     Ezt=Ez;
% % %     Hyt=Hy;
% % %     
% % %     Hy=Hy+dt/mu*(  ( Ez(2:X_num+1,:)-Ez(1:X_num,:) )./dx + ( Ex(:,1:Z_num)-Ex(:,2:Z_num+1) )./dz) ;
% % %     Ex(:,2:Z_num)=Ex(:,2:Z_num)+...
% % %         dt/epsilon*( ( Hy(:,1:Z_num-1)-Hy(:,2:Z_num) )/dz  );
% % %     Ez(2:X_num,:)=Ez(2:X_num,:)+...
% % %         dt/epsilon*( ( Hy(2:X_num,:)-Hy(1:X_num-1,:) )/dx );
% % %     
% % %     %% 边界条件
% % %     %   上下壁板
% % %     Ez(1,:)=0;
% % %     Ez(X_num+1,:)=0;
% % %     %   挡板
% % % %     Ex(1:fix(X_num/3),fix(Z_num/2)-2:fix(Z_num/2)+2)=0;
% % % %     Ex(fix(2*X_num/3):X_num,fix(Z_num/2)-2:fix(Z_num/2)+2)=0;
% % %     %
% % %     % 一阶Mur吸收条件 %吸收条件中真正起作用的其实是Ex，Ez的吸收
% % %     Ex(:,1)=Ext(:,2)+( (c*dt-dz)/(c*dt+dz)*( Ex(:,2)-Ext(:,1)) );
% % %     Ez(:,1)=Ezt(:,2)+( (c*dt-dz)/(c*dt+dz)*( Ez(:,2)-Ezt(:,1)) );
% % %     %     Hy(:,1)=Hyt(:,2)+( (c*dt-dz)/(c*dt+dz)*( Hy(:,2)-Hyt(:,1)) );
% % %     Ex(:,Z_num+1)=Ext(:,Z_num)+( (c*dt-dz)/(c*dt+dz)*( Ex(:,Z_num)-Ext(:,Z_num+1)) );
% % %     Ez(:,Z_num)=Ezt(:,Z_num-1)+( (c*dt-dz)/(c*dt+dz)*( Ez(:,Z_num-1)-Ezt(:,Z_num)) );
% % %     %     Hy(:,Z_num)=Hyt(:,Z_num-1)+( (c*dt-dz)/(c*dt+dz)*( Hy(:,Z_num-1)-Hyt(:,Z_num)) );
% % %     
% % %     %% 记录电压波形--屏蔽
% % %     %   TRM的设置
% % % %     TrmE_L1(:,:,t)=Ext(4:X_num-3,fix(Z_num/4));
% % % %     TrmE_R1(:,:,t)=Ext(4:X_num-3,fix(Z_num*3/4));
% % % %     Vref_Ltotal(t)=sum( Ext(4:X_num-3,fix(Z_num/4)) )*dx;
% % % %     Vref_Rtotal(t)=sum( Ext(4:X_num-3,fix(Z_num*3/4)) )*dx;
% % % %     Vref_BAN(t)=sum( Ext(4:X_num-3,fix(Z_num/2)) )*dx;
% % %     
% % %     %% 绘动图
% % %     PEx=Ex(1:X_num,1:Z_num);
% % %     PEz=Ez(1:X_num,1:Z_num);
% % %     Eabs=sqrt(PEx.^2+PEz.^2);
% % %     [xx,yy]=meshgrid(1:Z_num,1:X_num);
% % %     %     mesh(PEx)
% % %     mesh(xx,yy,PEx)
% % % %     view(0,90)
% % %     pause(0.00000001)
% % %     Tr_max=max(Ex,Ext);
% % % %     axis([])
% % % end
% % % 
% % % figure
% % % mesh(Tr_max)
