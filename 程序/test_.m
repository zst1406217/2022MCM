clc,clear
Japan=csvread('japan.csv',1,1);
alpha0=zeros(1,1050);
dx=22100/1047;
for i=1:1047
    alpha0(i+1)=atan((Japan(i+2)-Japan(i))/(2*dx));
end
xspan1=[0.01:0.01:5];
xspan2=[5.1:0.1:22100];
% x=[xspan1,xspan2];
x=[0.1:0.1:22100];
v=zeros(1,2710);
dt=zeros(1,2709);
t=zeros(1,2710);
P=300;
m=82.55;
CdA=0.21;
loss=0.02;
Crr=0.005;
Rho=1.22601;
k=0.5*CdA*Rho;
g=9.790;
v(1)=0.01;
for i=1:220999
    dl=(x(i+1)-x(i));
    v1=v(i);
    v(i+1)=sqrt((0.5*m*v1^2+P*dl/v1-(m*g*sin(alpha(x(i)))+Crr*m*g*cos(alpha(x(i)))+k*v1^2)*dl)/(0.5*m));
end
for i=1:220999
    dl=x(i+1)-x(i);
    dt(i)=dl/v(i);
    t(i+1)=t(i)+dt(i);
end
% plot(x,v,'-o');
plot(x,t,'-o');


function fun2=alpha(x)
global alpha0;
dx=22100/1047;
a1=(ceil(x/dx)-1)*dx;
a2=ceil(x/dx)*dx;
b1=alpha0(ceil(x/dx));
b2=alpha0(ceil(x/dx)+1);
fun2=(b2-b1)*(x-a1)/(a2-a1)+b1;
end