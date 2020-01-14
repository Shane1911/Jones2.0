%%�����غ������¾�������Qmax,delta_a,delta_r
clear all
close all
clc
global R_i alfa A Zn F K_oj K_ij f D gama omega_i m J dm Kn
% i=1;
% for Fa1=5000:5000:45000;
% Fa=2000;
% Fr=3000;
% M=150;%Mָ��y����
Fa=400;Fr=0;M=0;%Mָ��y����
omega_i=15000*pi/30;%rad/s
F=[Fa Fr M]';
% f=[0.519 0.519];%��������
% v=[0.3 0.3];    %���ɱ�
% dm=56.6;        %��Բֱ�� ��λ��mm
% D=13.5;         %����ֱ�� ��λ��mm
% alfa=5;        %��ʼ�Ӵ��� /���
% Zn=11;
%%���7008
f=[0.523 0.523];%��������
v=[0.3 0.3];    %���ɱ�
dm=53.9;        %��Բֱ�� ��λ��mm
D=7.3;         %����ֱ�� ��λ��mm
alfa=15;        %��ʼ�Ӵ��� /���
Zn=18;
%%�����ٶȶԱ����
% v=[0.3 0.3];%���ɱ�
% f=[0.53 0.56];%��������
% dm=54;%��λ��mm
% D=7.144;%��λ��mm
% alfa=15;%��ʼ�Ӵ���/���
% Zn=19;%���Ӹ���
%%���㾫ȷ�Ա�
% v=[0.3 0.3];%���ɱ�
% f=[0.523 0.523];%��������
% dm=125;%��λ��mm
% D=22.23;%��λ��mm
% alfa=40;%��ʼ�Ӵ���/���
% Zn=16;%���Ӹ���
%%wangB7014
% v=[0.3 0.3];%���ɱ�
% f=[0.523 0.523];%��������
% dm=89.5;%��λ��mm
% D=12.5;%��λ��mm
% alfa=25;%��ʼ�Ӵ���/���
% Zn=19;%���Ӹ���
%%��deta��J_a�������
alfa0=alfa*pi/180;
c=4.6568e-4;%�Ӵ�����ϵ��
fm=0.5*(f(1)+f(2));%ƽ����������
f1=@(x)(1+c/(2*fm-1)*(Fa/(Zn*D*D*sin(x)))^(2/3)-cos(alfa0)/cos(x));
alfa1=fzero(f1,alfa0+0.001);
%% �������
A_deta=[0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.25 1.67 2.5 5];
B_alfa=[1 0.9318 0.8964 0.8601 0.8225 0.7835 0.7427 0.6995 0.6529 0.6 0.4388 0.3088 0.185 0.0831];
C_J=[1/16 0.1707 0.211 0.2462 0.2782 0.3084 0.3374 0.3658 0.3945 0.4244 0.5044 0.606 0.724 0.8558];
alfa_deta=polyfit(B_alfa,A_deta,3);
alfa_J=polyfit(B_alfa,C_J,3);
x1=Fr*tan(alfa1)/Fa;
deta=polyval(alfa_deta,x1);
J_a=polyval(alfa_J,x1);
%��alfa1��������
f2=@(x)(1+c/(2*fm-1)*(1-1/2/(deta*cos(x)*cos(x)))*...
    (Fa/(Zn*D*D*J_a*sin(x)))^(2/3)-cos(alfa0)/cos(x));
alfa2=fzero(f2,alfa1+0.001);
%��������
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
%%����Qmax
Qmax=Fa/(Zn*J_a*sin(alfa2));%��λN
%����Kn
%%�նȼ���
int_K=[f dm D alfa v];
output_k=funK( int_K );
K2 =output_k(2);%��λ N*mm^1.5
K1 =output_k(1);%��λ N*mm^1.5
Kn=(1/((1/K1)^(2/3)+(1/K2^(2/3))))^1.5;%��λ N*m^1.5
%���㷨��������
delta_n=(Qmax/Kn)^(2/3);%��λ��mm
syms  xa yr
c1=sin(alfa2);
c2=cos(alfa2);
c3=tan(alfa2);
b1=delta_n;
b2=2*deta;
[xa,yr]=solve([xa*c1+yr*c2==b1,1+xa*c3/yr==b2],[xa,yr]);
delta_a=double(xa);% ��λ��mm
delta_r=double(yr);% ��λ��mm
output=[delta_a delta_r];
 %ϵͳ����:��Ȧ���أ�N N*m
R_i=(0.5*dm+(f(2)-0.5)*D*cosd(alfa));%�ڹ������������İ뾶/mm
A=(f(1)+f(2)-1)*D;%�������ľ���/m
 %��ֵ�趨:λ��/mm ��λ��/���
out=output;
x0=out(1);z0=out(2);theta=-0.01;
X0=[x0 z0 theta];
options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[X,fval,exitflag]=fsolve(@fun1,X0,options);
 if exitflag==1
     disp('��ֵ��Ч���������');
 else
     disp('�ó�ֵ�����޽�');
 end
  %�����⾲��ѧ��ֵ
   fan=0;
   A_1j=sind(alfa)+X(1)+R_i*X(3)*cos(fan);%��λ/mm
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
%ϵͳ����:��Ȧ���أ�N N*m
K_oj =K1;%��λ N*mm^1.5
K_ij =K2;%��λ N*mm^1.5
den=7850*1e-9;%�����ܶ� kg/mm^3
m=pi*den*D^3/6;%�������� kg
J=m*(D*1e-3)^2/10;%�����ת������ kg*m^2
gama=D/dm;
A=(f(1)+f(2)-1)*D;%�������ľ���/mm
 %��ֵ�趨:λ��/m ��λ�ƻ���
X3=output2;
Y3=output1;%��λ/mm
X1(1)=X3(1);X1(2)=X3(2);X1(3)=X3(3);
for kk=1:Zn
    X1(4+(kk-1)*4)=Y3(1);
    X1(5+(kk-1)*4)=Y3(2);
    X1(6+(kk-1)*4)=Y3(3);
    X1(7+(kk-1)*4)=Y3(4);
end
%��ʱ���⾲��ѧ����ʱ��
tic
 options =optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[Y,fval,exitflag]=fsolve(@fun2,X1,options);
if exitflag==1
    disp('�������');
else
    disp('�����ֵ��Ч');
end
toc
%��̬���ܲ������
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
   Q_1j(j)=K_oj*Y(6+(j-1)*4)^1.5;%��λ/N
   Q_2j(j)=K_ij*Y(7+(j-1)*4)^1.5;%��λ/N
   F_cj(j)=0.5*dm*(1e-3)*m*omega_cj^2;
   M_gj(j)=J*omega_Rj*omega_cj*sin_bj;
   F_1j=M_gj(j)*(1e3)*2/D;
   alfa_o(j)=acos(cos_1j)*180/pi;
   alfa_i(j)=acos(cos_2j)*180/pi;
   deta_o(j)=Y(6+(j-1)*4);
   deta_i(j)=Y(7+(j-1)*4);
end
%%������ȡ
% dateput(i,:)=[alfa_o(1) alfa_i(1) Q_1j(1) Q_2j(1) Y(1)];
% i=i+1;
% end
% save('date1.mat','dateput');
%%��ͼ
alfa_o(Zn+1)=alfa_o(1);
alfa_i(Zn+1)= alfa_i(1);
a=linspace(0,2*pi,Zn+1); %�趨�Ƕ�
figure(1)
b1=alfa_o; %�趨��Ӧ�Ƕȵİ뾶
b2=alfa_i; %�趨��Ӧ�Ƕȵİ뾶
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