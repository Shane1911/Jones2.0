%%联合载荷作用下静力计算Qmax,delta_a,delta_r
clear all
close all
clc
global R_i alfa A Zn F K_oj K_ij f D gama omega_i m J dm Kn
% i=1;
% for Fa1=5000:5000:45000;
% Fa=2000;
% Fr=3000;
% M=150;%M指向y方向
Fa=400;Fr=0;M=0;%M指向y方向
omega_i=15000*pi/30;%rad/s
F=[Fa Fr M]';
% f=[0.519 0.519];%沟道曲率
% v=[0.3 0.3];    %泊松比
% dm=56.6;        %节圆直径 单位：mm
% D=13.5;         %钢球直径 单位：mm
% alfa=5;        %初始接触角 /°度
% Zn=11;
%%轴承7008
f=[0.523 0.523];%沟道曲率
v=[0.3 0.3];    %泊松比
dm=53.9;        %节圆直径 单位：mm
D=7.3;         %钢球直径 单位：mm
alfa=15;        %初始接触角 /°度
Zn=18;
%%计算速度对比轴承
% v=[0.3 0.3];%泊松比
% f=[0.53 0.56];%沟道曲率
% dm=54;%单位：mm
% D=7.144;%单位：mm
% alfa=15;%初始接触角/°度
% Zn=19;%滚子个数
%%计算精确对比
% v=[0.3 0.3];%泊松比
% f=[0.523 0.523];%沟道曲率
% dm=125;%单位：mm
% D=22.23;%单位：mm
% alfa=40;%初始接触角/°度
% Zn=16;%滚子个数
%%wangB7014
% v=[0.3 0.3];%泊松比
% f=[0.523 0.523];%沟道曲率
% dm=89.5;%单位：mm
% D=12.5;%单位：mm
% alfa=25;%初始接触角/°度
% Zn=19;%滚子个数
%%对deta和J_a数据拟合
alfa0=alfa*pi/180;
c=4.6568e-4;%接触变形系数
fm=0.5*(f(1)+f(2));%平均沟道曲率
f1=@(x)(1+c/(2*fm-1)*(Fa/(Zn*D*D*sin(x)))^(2/3)-cos(alfa0)/cos(x));
alfa1=fzero(f1,alfa0+0.001);
%% 数据拟合
A_deta=[0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.25 1.67 2.5 5];
B_alfa=[1 0.9318 0.8964 0.8601 0.8225 0.7835 0.7427 0.6995 0.6529 0.6 0.4388 0.3088 0.185 0.0831];
C_J=[1/16 0.1707 0.211 0.2462 0.2782 0.3084 0.3374 0.3658 0.3945 0.4244 0.5044 0.606 0.724 0.8558];
alfa_deta=polyfit(B_alfa,A_deta,3);
alfa_J=polyfit(B_alfa,C_J,3);
x1=Fr*tan(alfa1)/Fa;
deta=polyval(alfa_deta,x1);
J_a=polyval(alfa_J,x1);
%求alfa1的收敛解
f2=@(x)(1+c/(2*fm-1)*(1-1/2/(deta*cos(x)*cos(x)))*...
    (Fa/(Zn*D*D*J_a*sin(x)))^(2/3)-cos(alfa0)/cos(x));
alfa2=fzero(f2,alfa1+0.001);
%迭代收敛
Rel=1;
while Rel>0.001
   x3=Fr*tan(alfa2)/Fa;
   deta=polyval(alfa_deta,x3);
   J_a=polyval(alfa_J,x3);
   f22=@(x)(1+c/(2*fm-1)*(1-1/2/(deta*cos(x)*cos(x)))*(Fa/(Zn*D*D*J_a*sin(x)))^(2/3)-cos(alfa0)/cos(x));
   alfa22=fzero(f22,alfa2+0.001);
   Rel=(abs(alfa22-alfa2))/alfa2;
   alfa2=alfa22;
end
%%计算Qmax
Qmax=Fa/(Zn*J_a*sin(alfa2));%单位N
%计算Kn
%%刚度计算
int_K=[f dm D alfa v];
output_k=funK( int_K );
K2 =output_k(2);%单位 N*mm^1.5
K1 =output_k(1);%单位 N*mm^1.5
Kn=(1/((1/K1)^(2/3)+(1/K2^(2/3))))^1.5;%单位 N*m^1.5
%计算法向趋近量
delta_n=(Qmax/Kn)^(2/3);%单位：mm
syms  xa yr
c1=sin(alfa2);
c2=cos(alfa2);
c3=tan(alfa2);
b1=delta_n;
b2=2*deta;
[xa,yr]=solve([xa*c1+yr*c2==b1,1+xa*c3/yr==b2],[xa,yr]);
delta_a=double(xa);% 单位：mm
delta_r=double(yr);% 单位：mm
output=[delta_a delta_r];
 %系统输入:内圈受载：N N*m
R_i=(0.5*dm+(f(2)-0.5)*D*cosd(alfa));%内滚道沟曲率中心半径/mm
A=(f(1)+f(2)-1)*D;%曲率中心距离/m
 %初值设定:位移/mm 角位移/°度
out=output;
x0=out(1);z0=out(2);theta=-0.01;
X0=[x0 z0 theta];
options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[X,fval,exitflag]=fsolve(@fun1,X0,options);
 if exitflag==1
     disp('初值有效，结果收敛');
 else
     disp('该初值附近无解');
 end
  %估算拟静力学初值
   fan=0;
   A_1j=sind(alfa)+X(1)+R_i*X(3)*cos(fan);%单位/mm
   A_2j=cosd(alfa)+X(2)*cos(fan);
   s=sqrt(A_1j^2+A_2j^2);
   Q=Kn*(A^1.5)*(s-1)^1.5;
   deta_oj=(Q/K1)^(2/3);
   deta_ij=(Q/K2)^(2/3);
   X1j=((f(1)-0.5)*D+deta_oj)*sind(alfa);
   X2j=(((f(1)-0.5)*D+deta_oj)^2-X1j^2)^0.5;
   Y1=[X1j X2j deta_oj deta_ij];
   output1=real(Y1);
   output2=real(X);
%    output1=[0.1869 0.4369 0.0077 0.0078];
%    output2=[0.0055 -0.0305 0.0039];
%系统输入:内圈受载：N N*m
K_oj =K1;%单位 N*mm^1.5
K_ij =K2;%单位 N*mm^1.5
den=7850*1e-9;%钢球密度 kg/mm^3
m=pi*den*D^3/6;%钢球质量 kg
J=m*(D*1e-3)^2/10;%钢球的转动惯量 kg*m^2
gama=D/dm;
A=(f(1)+f(2)-1)*D;%曲率中心距离/mm
 %初值设定:位移/m 角位移弧度
X3=output2;
Y3=output1;%单位/mm
X1(1)=X3(1);X1(2)=X3(2);X1(3)=X3(3);
for kk=1:Zn
    X1(4+(kk-1)*4)=Y3(1);
    X1(5+(kk-1)*4)=Y3(2);
    X1(6+(kk-1)*4)=Y3(3);
    X1(7+(kk-1)*4)=Y3(4);
end
%计时：拟静力学计算时间
tic
 options =optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[Y,fval,exitflag]=fsolve(@fun2,X1,options);
if exitflag==1
    disp('结果收敛');
else
    disp('所设初值无效');
end
toc
%动态性能参数输出
Q_1j=[];Q_2j=[];F_cj=[];M_gj=[];alfa_o=[]; alfa_i=[];
for j=1:Zn
   fan=(j-1)*2*pi/Zn;
   A_1j=A*sind(alfa)+Y(1)+R_i*Y(3)*cos(fan);
   A_2j=A*cosd(alfa)+Y(2)*cos(fan);
   cos_1j=Y(5+(j-1)*4)/((f(1)-0.5)*D+Y(6+(j-1)*4));
   sin_1j=Y(4+(j-1)*4)/((f(1)-0.5)*D+Y(6+(j-1)*4));
   cos_2j=(A_2j-Y(5+(j-1)*4))/((f(2)-0.5)*D+Y(7+(j-1)*4));
   sin_2j=(A_1j-Y(4+(j-1)*4))/((f(2)-0.5)*D+Y(7+(j-1)*4));
   tan_bj=sin_1j/(cos_1j+gama);
   cos_bj=1/sqrt(1+tan_bj^2);
   sin_bj=tan_bj*cos_bj;
   omega_cj=omega_i*(1-gama*cos_2j)/(1+cos_2j*cos_1j+sin_1j*sin_2j);
   omega_Rj=-omega_i/(((cos_1j+tan_bj*sin_1j)/(1+gama*cos_1j)+(cos_2j+tan_bj*sin_2j)/(1-gama*cos_2j))*gama*cos_bj);
   Q_1j(j)=K_oj*Y(6+(j-1)*4)^1.5;%单位/N
   Q_2j(j)=K_ij*Y(7+(j-1)*4)^1.5;%单位/N
   F_cj(j)=0.5*dm*(1e-3)*m*omega_cj^2;
   M_gj(j)=J*omega_Rj*omega_cj*sin_bj;
   F_1j=M_gj(j)*(1e3)*2/D;
   alfa_o(j)=acos(cos_1j)*180/pi;
   alfa_i(j)=acos(cos_2j)*180/pi;
   deta_o(j)=Y(6+(j-1)*4);
   deta_i(j)=Y(7+(j-1)*4);
end
%%数据提取
% dateput(i,:)=[alfa_o(1) alfa_i(1) Q_1j(1) Q_2j(1) Y(1)];
% i=i+1;
% end
% save('date1.mat','dateput');
%%作图
alfa_o(Zn+1)=alfa_o(1);
alfa_i(Zn+1)= alfa_i(1);
a=linspace(0,2*pi,Zn+1); %设定角度
figure(1)
b1=alfa_o; %设定对应角度的半径
b2=alfa_i; %设定对应角度的半径
polar(a,b2,'-k^')
hold on
polar(a,b1,'-r^')
hold on
grid on
view([-90 90])
legend('contact angle \alpha_i','contact angle \alpha_o') 
figure(2)
subplot(3,2,1)
plot(1:Zn,alfa_o(1:Zn),'-r^');
hold on
grid on
plot(1:Zn,alfa_i(1:Zn),'-k^');
legend(' \alpha_o','\alpha_i') 
subplot(3,2,2)
plot(1:Zn,Q_1j,'-r^');
hold on
grid on
plot(1:Zn,Q_2j,'-k^');
legend('Contact force Q_o','Contact force Q_i') 
subplot(3,2,3)
plot(1:Zn,deta_o,'-r^');
hold on
grid on
plot(1:Zn,deta_i,'-k^');
legend('outer \delta_o','inner \delta_i') 
subplot(3,2,4)
plot(1:Zn,F_cj,'-r^');
hold on
grid on
legend('Centrifugal force F_c') 
subplot(3,2,5)
plot(1:Zn,M_gj,'-r^');
hold on
grid on
legend('Gyroscopic moment M_{gy}') 
figure(3)
plot(1:Zn,Q_1j,'-r^');
hold on
plot(1:Zn,Q_2j,'-k^');
legend('Contact force Q_o','Contact force Q_i') 