# %%
#  python -i japan.py | tee out.txt 
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from scipy import  interpolate 
# %%
k_=20
Japan=pd.read_csv('bel-women-h.csv')
del Japan['Unnamed: 0']
Japan=np.array(Japan).reshape(len(Japan))
# %%
len(Japan)
s_=44300
dx=s_/1461
m=82.55
CdA=0.21
loss=0.02
Crr=0.005
Rho=1.22601
k=0.5*CdA*Rho
g=9.790    
W=10860#2.079e4
CP=299.6#422.8
dx_P=s_/(100*k_)
alpha0=np.zeros(1090)
y__0=np.zeros(1090)
# dx=s_/1048
alpha0[1:1087]=(Japan[:-2]-Japan[2:])/2/dx
alpha0[0]=alpha0[1]
y__0[1:1087]=(Japan[2:]+Japan[:-2]-2*Japan[1:-1])/((dx)**2)
y__0[0]=y__0[1]
start=np.full(k_*100,299.6)
lb=np.ones(k_*100)
ub=np.full(k_*100,1706)
T=[]
X=[]
# %%
P_=interpolate.interp1d(np.linspace(0,s_,100*k_),start,kind='cubic')
alpha=interpolate.interp1d(np.linspace(0,s_,1086),alpha0[:1086],kind='cubic')
y__=interpolate.interp1d(np.linspace(0,s_,1086),y__0[:1086],kind='cubic')
# %%
def constrains(P0):
    AAW=np.zeros(k_*100)
    AAW[0]=W
    g=np.zeros(200*k_)
    for i in range(100*k_-1):
        g[2*i]=W-AAW[i]
        g[2*i+1]=AAW[i]
        if(g[2*i]<0 or g[2*i+1]<0): 
            return 0
        l=np.searchsorted(X,i*dx_P,'right')
        r=np.searchsorted(X,i*dx_P,'left')
        P_average=P0[i]
        if P_average>=CP :
            # AAW[i+1]=AAW[i]-(P_average-CP)*(T[0][r]-T[0][l])
            AAW[i+1]=AAW[i]-(P_average-CP)*(np.sum(T[l:r+1]))
        elif P_average<CP:
            taoW=546*np.exp(-0.01*(CP-P_average))+316
            AAW[i+1]=W-(W-AAW[i])*np.exp(-np.sum(T[l:r+1])/taoW)
            # AAW[i+1]=W-(W-AAW[i])*np.exp(T[0][l]-T[0][r]/taoW)
    g[(200*k_-2)]=W-AAW[100*k_-1]
    g[(200*k_-1)]=AAW[100*k_-1]
    if(g[200*k_-2]<0 or g[200*k_-1]<0): 
        return 0
    # h=0
    return 1
res=[]
min=1e4
def totalObject(P0):
    pre=0
    global T 
    global X
    T=[]
    X=[]
    P0=P0*(1-loss)
    P_=interpolate.interp1d(np.linspace(0,s_,100*k_),P0,kind='cubic')
    def subObject(l,r,n,P0):
        def func(x,s):
            alpha_x=alpha(x+l)
            return (-P_(x+l)/m*s**4*(np.cos(alpha_x))**2+k/(m*np.cos(alpha_x))*s\
                +s*(g*s**2-y__(x+l))*(Crr+np.tan(alpha_x)*np.cos(2*alpha_x)))
        y0=[0.4]
        y3 = integrate.solve_ivp(func,t_span=(0,r-l), y0=y0, method='RK45',\
             t_eval=np.linspace(0,r-l,n)) 
        sub_sum=np.sum(np.dot(y3.y[0][1:]+y3.y[0][:-1],y3.t[1:]-y3.t[:-1]))/2
        if pre:
            print('l %d r %d n %d N %d sub sum %f'%(l,r,n,y3.t.shape[0],sub_sum))
        return (sub_sum,(y3.y[0]).tolist(),(y3.t).tolist())
    sum=0
    for i in range(30):
        subsum,T_,X_=subObject(1000*i,1000*(i+1),50,P0)
        T+=T_
        X+=X_
        if subsum<0:
            sum=4e4+np.random.rand()*100

            if pre:
                print(sum)
            return sum
        else :
            sum+=subsum
    print('sum',sum)
    T=np.array(T)
    X=np.array(X)
    if(constrains(P0)==0):
        sum=40000
        print(sum)
    if(sum<2160):
        init=pd.DataFrame(P0)
        init.to_csv('bel-women-init-2160.csv')
    global res 
    global min
    if sum<min:
        min=sum 
        res=P0
    else:
        print('min',min)
    return sum

# %%
np.random.seed(42)
start=pd.read_csv('bel-women-init-2178.csv')['0'].to_numpy()
res = minimize(totalObject,start,method='Nelder-Mead',bounds=np.vstack([lb,ub]).T)
