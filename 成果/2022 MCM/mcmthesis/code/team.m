clc,clear
CP_High=422.8;W_High=20790;CP_Low=412;W_Low=18000;P=300;m=82.55;CdA=0.21;loss=0.02;
Crr=0.005;Rho=1.22601;k=0.5*CdA*Rho;g=9.790;
rate=[0.975 0.616 0.497 0.443 0.426 0.433];
t_High=20;t_Low=10;start=zeros(1,9);
for i=1:6
    start(i)=k*rate(i)*16^3+m*g*Crr*16;
end
start(7)=k*(16-1.17)^3+m*g*Crr*(16-1.17);
start(8)=16;
start(9)=15;
x = fsolve(@(x) team_fun(x,CP_High,CP_Low,k,rate,m,Crr,g),start) 
t_Low=x(9);
V_team=x(8);
V_straight=(V_team*(3*(20+t_Low)-10)+(V_team-1.17)*10)/(3*(20+t_Low));
t_High=[0,20,30,20+x(9),40+x(9),40+2*x(9),60+2*x(9),60+3*x(9)];
AAW_High=zeros(1,8);
AAW_High(1)=W_High;
for i=2:8
    AAW_High(i)=AAW_High(i-1)+(CP_High-x(i-1))*(t_High(i)-t_High(i-1));
end
plot(t_High,AAW_High,'-*')