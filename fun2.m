function [ Y ] = fun2( X );
global R_i alfa A Zn F K_oj K_ij f D gama omega_i m J dm  
Y=zeros(4*Zn+3,1);
for j=1:Zn
    e=0;
    fan=(j-1)*2*pi/Zn;
    A_1j=A*sind(alfa)+X(1)+R_i*X(3)*cos(fan);
    A_2j=A*cosd(alfa)+X(2)*cos(fan)+e;
    cos_1j=X(5+(j-1)*4)/((f(1)-0.5)*D+X(6+(j-1)*4)-e);
    sin_1j=X(4+(j-1)*4)/((f(1)-0.5)*D+X(6+(j-1)*4)-e);
    cos_2j=(A_2j-X(5+(j-1)*4))/((f(2)-0.5)*D+X(7+(j-1)*4)-e);
    sin_2j=(A_1j-X(4+(j-1)*4))/((f(2)-0.5)*D+X(7+(j-1)*4)-e);
    tan_bj=sin_1j/(cos_1j+gama);
    cos_bj=cos(atan(tan_bj));
    sin_bj=tan_bj*cos_bj;
    omega_cj=omega_i*(1-gama*cos_2j)/(1+cos_2j*cos_1j+sin_1j*sin_2j);
    omega_Rj=-omega_i/(((cos_1j+tan_bj*sin_1j)/(1+gama*cos_1j)+...
        (cos_2j+tan_bj*sin_2j)/(1-gama*cos_2j))*gama*cos_bj);
    Q_1j=K_oj*X(6+(j-1)*4)^1.5;%单位/N
    Q_2j=K_ij*X(7+(j-1)*4)^1.5;%单位/N
    F_cj=0.5*dm*(1e-3)*m*omega_cj^2;
    M_gj=J*omega_Rj*omega_cj*sin_bj;
    F_1j=M_gj*(1e3)*2/D;
    Y(4*j-3,1)=-Q_1j*sin_1j+Q_2j*sin_2j-F_1j*cos_1j;
    Y(4*j-2,1)=-Q_1j*cos_1j+Q_2j*cos_2j+F_1j*sin_1j+F_cj;
    Y(4*j-1,1)=(A_1j-X(4+(j-1)*4))^2+(A_2j-X(5+(j-1)*4))^2-((f(2)-0.5)*D+X(7+(j-1)*4)-e)^2;
    Y(4*j,1)=X(4+(j-1)*4)^2+X(5+(j-1)*4)^2-((f(1)-0.5)*D+X(6+(j-1)*4)-e)^2;
    Y(4*Zn+1,1)=Y(4*Zn+1,1)+Q_2j*sin_2j;
    Y(4*Zn+2,1)=Y(4*Zn+2,1)+Q_2j*cos_2j*cos(fan);
    Y(4*Zn+3,1)=Y(4*Zn+3,1)+Q_2j*R_i*1e-3*sin_2j*cos(fan);
end
Y(4*Zn+1,1)=Y(4*Zn+1,1)-F(1);
Y(4*Zn+2,1)=Y(4*Zn+2,1)-F(2);
Y(4*Zn+3,1)=Y(4*Zn+3,1)-F(3);
end

