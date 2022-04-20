clc,clear
a_t=0.99;
warning on;
tf=3;t_current=97;
Markov_length=100000;
Japan=csvread('japan.csv',1,1);
alpha0=zeros(1,1050);
y__0=zeros(1,1050);
dx=22100/1047;
m=82.55;
CdA=0.21;
loss=0.02;
Crr=0.005;
Rho=1.22601;
k=0.5*CdA*Rho;
g=9.790;
xspan=[0.01 22080];
y0=10;
for i=1:1047
    alpha0(i+1)=atan((Japan(i+2)-Japan(i))/(2*dx));
    y__0(i+1)=(Japan(i+2)+Japan(i)-2*Japan(i+1))/((dx)^2);
end
sol_current=repelem(300,7);
sol_new=sol_current;
sol_best=sol_current;
T_current=inf;T_best=inf;
while t_current>tf
    for r=1:Markov_length
        clear v s x t t0 x0 g0;
%         if(rand<0.5)
%             ind=floor(rand*99+1);
%             while sol_new(ind)>=590
%                 ind=floor(rand*100);
%             end
%             sol_new(ind)=sol_new(ind)+3;
%         else
%             ind=floor(rand*99+1);
%             while sol_new(ind)<=9
%                 ind=floor(rand*100);
%             end
%             sol_new(ind)=sol_new(ind)-3;
%         end
        sol_new=zeros(1,7);
        for i=1:3
            sol_new(i*2-1)=unidrnd(300)+300;
        end
        sol_new(2)=unidrnd(300);
        sol_new(4)=unidrnd(300);
        sol_new(6:7)=300;
        flag=1;
        x=[0.1:0.1:22100];
        v=zeros(1,221000);
        dt=zeros(1,221000);
        t=zeros(1,221000);
        v(1)=0.01;
        for i=1:220999
            dl=(x(i+1)-x(i))/cos(alpha(x(i),dx,alpha0));
            v1=v(i);
            if (0.5*m*v1^2+P(x(i),sol_new)*dl/v1-(m*g*sin(alpha(x(i),dx,alpha0))+Crr*m*g*cos(alpha(x(i),dx,alpha0))+k*v1^2)*dl)/(0.5*m)<0
                flag=0;
                break;
            end
            v(i+1)=sqrt((0.5*m*v1^2+P(x(i),sol_new)*dl/v1-(m*g*sin(alpha(x(i),dx,alpha0))+Crr*m*g*cos(alpha(x(i),dx,alpha0))+k*v1^2)*dl)/(0.5*m));
        end
        if flag==1
            for i=1:220999
                dl=x(i+1)-x(i);
                dt(i)=dl/v(i);
                t(i+1)=t(i)+dt(i);
            end
            %         [x,s]=ode45(@(x,s) -P(x,sol_new)/m*s^4*(cos(alpha(x,dx,alpha0)))^2+k/(m*cos(alpha(x,dx,alpha0)))*s+s*(g*s^2-y__(x,dx,y__0))*(Crr+tan(alpha(x,dx,alpha0))*cos(2*(alpha(x,dx,alpha0)))),xspan,y0);
            %         t(1)=0;
            size_x=size(x);
            %         if isnan(s(size_x(1)))
            %             flag=0;
            %         end
            %         for i=1:size_x(1)-1
            %             t(i+1)=t(i)+(x(i+1)-x(i))*s(i);
            %         end
            t0=zeros(1,7);
            x0=dx*[0 82 366 617 793 845 1047];
            t0(1)=0;
            t0(7)=t(size_x(2));
            for i=2:6
                for j=1:size_x(2)-1
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
            W=10860;
            AAW=zeros(1,7);
            AAW(1)=W;
            CP=299.6;
            g0=zeros(1,14);
            for i=1:6
                g0(2*i-1)=AAW(i)-W;
                g0(2*i)=-AAW(i);
                P_average=(sol_new(i)+sol_new(i+1))/2;
                if P_average>=CP
                    AAW(i+1)=AAW(i)-(P_average-CP)*(t0(i+1)-t0(i));
                elseif P_average<CP
                    AAW(i+1)=AAW(i)+(CP-P_average)*(t0(i+1)-t0(i));
                end
            end
            g0(13)=AAW(7)-W;
            g0(14)=-AAW(7);
            for i=1:14
                if g0(i)>0
                    flag=0;
                    break;
                end
            end
            T_new=t(size_x(2));
            if flag==1
                if T_new<T_current
                    T_current=T_new;
                    sol_current=sol_new;
                    if T_new<T_best
                        T_best=T_new
                        sol_best=sol_new;
                    end
                else
                    if rand<exp(-(T_new-T_current)/t_current)
                        T_current=T_new;
                        sol_current=sol_new;
                    else
                        sol_new=sol_current;
                    end
                end
            end
        end
    end
    t_current=t_current*a_t
end
disp(sol_best)
disp(T_best)



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
