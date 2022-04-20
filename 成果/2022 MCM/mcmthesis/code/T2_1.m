clc,clear
global Japan alpha0 y__0;
warning off;
Japan=csvread('japan.csv',1,1);
alpha0=zeros(1,1050);
y__0=zeros(1,1050);
dx=22100/1047;
for i=1:1047
    alpha0(i+1)=atan((Japan(i+2)-Japan(i))/(2*dx));
    y__0(i+1)=(Japan(i+2)+Japan(i)-2*Japan(i+1))/((dx)^2);
end
start=repelem(299.6,5);
lb=zeros(1,5);
ub=repelem(600,5);
options = optimoptions('fmincon','MaxIterations',10000,'MaxFunctionEvaluations',...
500000,'StepTolerance',1.00000e-5,'ConstraintTolerance',1.00000e-7,'Algorithm','sqp');
[P,fval]=fmincon('T2_1_object',start,[],[],[],[],lb,ub,'T2_1_constrains',options)