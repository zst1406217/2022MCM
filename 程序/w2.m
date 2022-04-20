%两段路，无坡度，前一半路程逆风，后一半顺风，均匀速，考虑恢复
clc;clear;
length=2000;
P_a=10000;
Crr=0.0042; mg=800;
CP=300;
v_1=5;
v_2=-5;
k=0.12;
b1=0.089;
b2=0.0084;
syms x;
T_min=200000;
% T_f=zeros(8,1);
% P=zeros(8,2);
turns=1;

for P1=300:0.5:320 %p1 逆风 
     W_ex=0;
     W=P_a;
     v_t=double(solve(k*(x+5)^3+Crr*mg*x+b1*x+b2*x^2-P1, 'Real', true));
     t1=length/2/v_t;
      if (P1-CP)*t1>P_a
          break;
      end
      if P1>CP
        W=P_a-(P1-CP)*t1;
        W_ex=W_ex+(P1-CP)*t1;
      else 
         W=P_a-W_ex*exp(-(P1-CP)/P_a);
      end
      T_min2=200000;
      for P2=300:0.1:320 %精度加细
          v_t2=double(solve(k*(x-5)^3+Crr*mg*x-P2, 'Real', true));
          t2=length/2/v_t2;
          if(P2-CP)*t2>W
              break;
          end
          T=t1+t2;
          if T<T_min2
              T_min2=T;
              P(turns,2)=P2;
          end
      end
      T_f(turns)=T_min2;
      P(turns,1)=P1;
      temp = [T_min2 P1 P2]
      turns=turns+1;
       if T_min2<T_min
              T_min=T_min2;
       end
end
plot(P(:,1),T_f)