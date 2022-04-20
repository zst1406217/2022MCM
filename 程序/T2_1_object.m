function f = T2_1_object(P0)
global alpha0 y__0;
P0=[P0 299.6 299.6];
s=22100;
dx=s/1047;
m=82.55;
CdA=0.21;
loss=0.02;
Crr=0.005;
Rho=1.22601;
k=0.5*CdA*Rho;
g=9.790;
P0=P0*(1-loss);
xspan=[0.01 22080];
y0=10;
[x,s]=ode23s(@(x,s) -P(x,P0)/m*s^4*(cos(alpha(x,dx,alpha0)))^2+k/(m*cos(alpha(x,dx,alpha0)))*s+s*(g*s^2-y__(x,dx,y__0))*(Crr+tan(alpha(x,dx,alpha0))*cos(2*(alpha(x,dx,alpha0)))),xspan,y0);
t(1)=0;
size_x=size(x);
for i=1:size_x(1)-1
    t(i+1)=t(i)+(x(i+1)-x(i))*s(i);
end
x0=dx*[0 82 366 617 793 845 1047];
t0=zeros(1,7);
t0(7)=t(size_x(1));
for i=2:6
    for j=1:size_x(1)-1
        if x(j)<=x0(i) && x(j+1)>x0(i)
            j0=j;
            break;
        end
    end
    a1_=x(j0);
    a2_=x(j0+1);
    b1_=t(j0);
    b2_=t(j0+1);
    t0(i)=(b2_-b1_)*(x0(i)-a1_)/(a2_-a1_)+b1_;
end
output=[x0;t0];
save output;
sum=t(size_x(1));
f=sum
end

function fun1=P(x,P0)
dx=22100/1047;
x0=dx*[0 82 366 617 793 845 1047];
for i=1:7
    if x0(i)<=x && x0(i+1)>x
        i0=i;
        break;
    end
end
a1=x0(i);
a2=x0(i+1);
b1=P0(i);
b2=P0(i+1);
fun1=(b2-b1)*(x-a1)/(a2-a1)+b1;
end

function fun2=alpha(x,dx,alpha0)
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


