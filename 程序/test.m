clc,clear
Japan=csvread('japan.csv',1,1);
global alpha0 dx;
% P=488;
s=22100;
dx=s/1047;
m=82.55;
CdA=0.21;
loss=0.02;
Crr=0.005;
Rho=1.22601;
k=0.5*CdA*Rho;
s=zeros(1,1050);
s_=zeros(1,1050);
% dt=zeros(1,1050);
g=9.790;
% P=P*(1-loss);
alpha0=zeros(1,1050);
y__0=zeros(1,1050);
for i=1:1047
    alpha0(i+1)=atan((Japan(i+2)-Japan(i))/(2*dx));
    y__0(i+1)=(Japan(i+2)+Japan(i)-2*Japan(i+1))/((dx)^2);
end
xspan1=[0.01:0.01:5];
xspan2=[6:10:22080];
xspan=[xspan1,xspan2];
y0=10;
[x,s]=ode23(@(x,s) -P(x)/m*s^4*(cos(alpha(x)))^2+k/(m*cos(alpha(x)))*s+s*(g*s^2-y__(x,dx,y__0))*(Crr+tan(alpha(x))*cos(2*(alpha(x)))),xspan,y0);
t(1)=0;
size_x=size(x);
for i=1:size_x(1)-1
    t(i+1)=t(i)+(x(i+1)-x(i))*s(i);
end
t0=zeros(1,100);
x0=zeros(1,100);
t0(1)=0;
dx_t0=22080/99;
size_x=size(x);
for i=1:99
    for j=1:size_x(1)-1
        if (x(j)<=i*dx_t0)&&(x(j+1)>=i*dx_t0)
            j0=j;
            break;
        end
    end
    
    a1_=x(j0);
    a2_=x(j0+1);
    b1_=t(j0);
    b2_=t(j0+1);
    t0(i+1)=(b2_-b1_)*(i*dx_t0-a1_)/(a2_-a1_)+b1_;
    x0(i+1)=i*dx_t0;
end
output=[x0;t0];
save output;
t=t(1:size(x));
size=size(x);
sum=t(size(1));
plot(x,t,'-o')

function fun1=P(x)
if x<=800
    fun1=299;
else
    fun1=299;
end
end

function fun2=alpha(x)
global dx alpha0;
a1=(ceil(x/dx)-1)*dx;
a2=ceil(x/dx)*dx;
b1=alpha0(ceil(x/dx));
b2=alpha0(ceil(x/dx)+1);
fun2=(b2-b1)*(x-a1)/(a2-a1)+b1;
end

function fun3=y__(x,dx,y__0)
a1=(ceil(x/dx)-1)*dx;
a2=ceil(x/dx)*dx;
b1=y__0(ceil(x/dx));
b2=y__0(ceil(x/dx)+1);
fun3=(b2-b1)*(x-a1)/(a2-a1)+b1;
end