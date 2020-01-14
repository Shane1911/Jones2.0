%%数据整理
clear all
close all
clc
P1=load('date1.mat','dateput');
date=P1.dateput;
%%接触角变化
X=[0 5 10 15 20 25 30 35 40 45];
Y1=[59 57.5 55 53 49 48 47.5 47 46.5 47];
Y2=[0 10 19 26 32 34 36 37 37.5 38];
Qo=[2400 2450 2500 2700 2850 3050 3400 3800 4200 4400];
Qi=[0 450 700 1000 1500 1900 2250 2600 3100 3500];
delta=[-0.25 -0.12 -0.045 0 0.02 0.035 0.045 0.05 0.055 0.06];
%计算数据提取
Y3=[];Y4=[];Q3=[];Q4=[];delta_a=[];
Y3(1)=60;Q3(1)=2400;
Y4(1)=0;Q4(1)=0;
delta_a(1)=-0.25;
n=length(date(:,1));
for kk=1:n
    Y3(kk+1)=date(kk,2);
    Y4(kk+1)=date(kk,1);
    Q3(kk+1)=date(kk,3);
    Q4(kk+1)=date(kk,4);
    delta_a(kk+1)=date(kk,5);
end
P1=polyfit(X,Y1,3);
P2=polyfit(X,Y2,3);
P3=polyfit(X,Qo,4);
P4=polyfit(X,Qi,4);
P5=polyfit(X,delta,4);
x1=0:5:45;
y1=polyval(P1,x1);
y2=polyval(P2,x1);
Q1=polyval(P3,x1);
Q2=polyval(P4,x1);
Delta=polyval(P5,x1);
%可决系数
R_alfa_o=sqrt((sum(y1.*Y3))^2/(sum(y1.^2)*sum(Y3.^2)))
R_alfa_i=sqrt((sum(y2.*Y4))^2/(sum(y2.^2)*sum(Y4.^2)))
R_delta_a=sqrt((sum(Delta.*delta_a))^2/(sum(Delta.^2)*sum(delta_a.^2)))
figure(1)
h1=plot(x1,y1,'k-','LineWidth',1.5);
hold on
h3=plot(x1,y2,'r-','LineWidth',1.5);
hold on
h2=plot(X,Y3,'k^','MarkerSize',9);
hold on
h4=plot(X,Y4,'r^','MarkerSize',9);
axis([0 47 0 65])
xlabel('轴向载荷 F_a,kN')
ylabel('接触角 \alpha_i and \alpha_o')
legend([h1 h2],'Harris 计算结果','本发明结果')
grid on
hold on
figure(2)
h1=plot(x1,Q1,'k-','LineWidth',1.5);
hold on
h3=plot(x1,Q2,'r-','LineWidth',1.5);
hold on
h2=plot(X,Q3,'k^','MarkerSize',9);
hold on
h4=plot(X,Q4,'r^','MarkerSize',9);
xlabel('Axial load F_a,kN')
ylabel('Contact angle \alpha_i and \alpha_o')
legend([h1 h2],'Harris 计算结果','本发明结果')
grid on
hold on
figure(3)
h1=plot(x1,Delta,'k-','LineWidth',1.5);
hold on
h2=plot(x1,delta_a,'r^','LineWidth',1.5);
xlabel('轴向载荷 F_a,kN')
ylabel('轴向位移 \delta_a /mm')
legend([h1 h2],'Harris 计算结果','本发明结果')
grid on
hold on