%%��Բ��������
clear all
close all
clc
f=[0.53 0.56];%��������
dm=54;%��λ��mm
D=7.144;%��λ��mm
alfa=15;%��ʼ�Ӵ���/���
v=[0.3 0.3];
alfa0=15;
R_delta1=0.5*(dm/cosd(alfa0)+D);
R=R_delta1*D/(R_delta1+D);
E12b=1e5*[2.08 2.08];%��λ��N/mm^2
E=2/((1-v(1)^2)/E12b(1)+(1-v(2)^2)/E12b(2))
R_delta1=0.5*D*(dm+D*cosd(alfa0))/dm;
R_delta2=0.5*D*(dm-D*cosd(alfa0))/dm;
R_yita1=f(1)*D/(2*f(1)-1);
R_yita2=f(2)*D/(2*f(2)-1);
deta1=1.0003+0.5968/(R_yita1/R_delta1)
deta2=1.0003+0.5968/(R_yita2/R_delta2)
tao1=1.5277-0.6023*log(R_delta1/R_yita1);
tao2=1.5277-0.6023*log(R_delta2/R_yita2);
R1=(R_yita1*R_delta1)/(R_delta1+R_yita1)%��λ��mm
R2=(R_yita2*R_delta2)/(R_delta2+R_yita2)%��λ��mm
k1=1.0339*(R_yita1/R_delta1)^0.636
k2=1.0339*(R_yita2/R_delta2)^0.636
Ko=pi*k1*E*sqrt(R1*deta1/(4.5*(tao1)^3))%N*mm^1.5
Ki=pi*k2*E*sqrt(R2*deta1/(4.5*(tao2)^3))%N*mm^1.5
%a=(6*k1^2*deta1*Q_j*R1/(pi*E))^(1/3)