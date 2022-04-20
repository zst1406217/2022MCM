clc,clear
rand('seed',sum(100*clock));
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
ttt=[340 260 331 269 328 300 300];
sol_1=[ttt(1):(ttt(2)-ttt(1))/7:ttt(2)];
sol_2=[ttt(2)+(ttt(3)-ttt(2))/8:(ttt(3)-ttt(2))/8:ttt(3)];
sol_3=[ttt(3)+(ttt(4)-ttt(3))/45:(ttt(4)-ttt(3))/45:ttt(4)];
sol_4=[ttt(4)+(ttt(5)-ttt(4))/18:(ttt(5)-ttt(4))/18:ttt(5)];
sol_5=[ttt(5)+(ttt(6)-ttt(5))/5:(ttt(6)-ttt(5))/5:300];
sol_6=repelem(300,16);
sol_current=[sol_1 sol_2 sol_3 sol_4 sol_5 sol_6];
sol_current=floor(sol_current)
sol_new=sol_current;
sol_best=sol_current;
T_current=inf;T_best=inf;
while t_current>tf
    for r=1:Markov_length
        clear v s x t t0 x0 g0;
        sol_new=zeros(1,100);
        for i=1:100
            sol_new(i)=unidrnd(600);
        end
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
            t0=zeros(1,100);
            x0=zeros(1,100);
            t0(1)=0;
            dx_t0=22080/99;
            
            for i=1:99
                for j=1:size_x(2)-1
                    if (x(j)<=i*dx_t0)&&(x(j+1)>i*dx_t0)
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
            W=10860;
            AAW=zeros(1,100);
            AAW(1)=W;
            CP=299.6;
            g0=zeros(1,200);
            for i=1:99
                g0(2*i-1)=AAW(i)-W;
                g0(2*i)=-AAW(i);
                P_average=(sol_new(i)+sol_new(i+1))/2;
                if P_average>=CP
                    AAW(i+1)=AAW(i)-(P_average-CP)*(t0(i+1)-t0(i));
                elseif P_average<CP
                    taoW=546*exp(-0.01*(CP-P_average))+316;
                    AAW(i+1)=W-(W-AAW(i))*exp(-(t0(i+1)-t0(i))/taoW);
                end
            end
            g0(199)=AAW(100)-W;
            g0(200)=-AAW(100);
            for i=1:200
                if g0(i)>0
                    flag=0
                    break;
                end
            end
            T_new=t(size_x(2))
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
dx_P=22100/99;
a1=(ceil(x/dx_P)-1)*dx_P;
a2=ceil(x/dx_P)*dx_P;
b1=P0(ceil(x/dx_P));
b2=P0(ceil(x/dx_P)+1);
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