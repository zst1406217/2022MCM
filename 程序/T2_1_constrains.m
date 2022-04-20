function [g,h] = T2_1_constrains(P0)
load('output.mat');
P0=[P0 299.6 299.6]
t0=output(2,:);
W=10860;
AAW=zeros(1,7);
AAW(1)=W;
CP=299.6;
g=zeros(1,14);
for i=1:6
   g(2*i-1)=AAW(i)-W;
   g(2*i)=-AAW(i);
   P_average=(P0(i)+P0(i+1))/2;
   if P_average>=CP
       AAW(i+1)=AAW(i)-(P_average-CP)*(t0(i+1)-t0(i));
   elseif P_average<CP
%        taoW=546*exp(-0.01*(CP-P_average))+316;
%        AAW(i+1)=W-(W-AAW(i))*exp(-(t0(i+1)-t0(i))/taoW);
       AAW(i+1)=AAW(i)+(CP-P_average)*(t0(i+1)-t0(i));
   end
end
g(13)=AAW(7)-W;
g(14)=-AAW(7);
h=0;
save g;
end

