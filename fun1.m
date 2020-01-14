function [ Y ] = fun1( X );
global R_i alfa Kn A Zn F
X=X/A;%无量纲化数据
F_L=zeros(3,1);
%计算求和方程
for j=1:Zn
    %相关参数计算
    e=0;
    fan=(j-1)*2*pi/Zn;
    x1=X(1);x2=X(2);x3=X(3);
    A_1j=sind(alfa)+X(1)+R_i*X(3)*cos(fan);%单位/mm
    A_2j=cosd(alfa)+X(2)*cos(fan)+e;
    s=sqrt(A_1j^2+A_2j^2);
    sin_alfa1=A_1j/s;
    cos_alfa1=A_2j/s;
    Q=Kn*(A^1.5)*(s-1)^1.5;
    %方程负载区求和
    F_L(1,1)=F_L(1,1)+Q*sin_alfa1;
    F_L(2,1)=F_L(2,1)+Q*cos(fan)*cos_alfa1;
    F_L(3,1)=F_L(3,1)+R_i*1e-3*Q*cos(fan)*sin_alfa1;
end
%目标函数
Y(1)=F_L(1,1)-F(1);
Y(2)=F_L(2,1)-F(2);
Y(3)=F_L(3,1)-F(3);
end
