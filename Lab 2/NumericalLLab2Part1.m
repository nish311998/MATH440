%C.2.1 Numerical Solution of IVP

%C2.1.1 Accuracy of RK Method
Nval=[10;20;40;80;160]; %Our different constant stepsizes
%y'=v=f(t,y,v)
%v'=(y^2-1)v-y=g(t,y,v)

e=zeros(1,5); %error vector

%Finding Nmax to use for error

%These are our intital values
Nmaxt=0; Nmaxtend=1; NmaxT=Nmaxt;
Nmax=320;
Nmaxh=Nmaxtend/Nmax;
Nmaxy0=1; Nmaxy(1)=Nmaxy0;
Nmaxdy0=0; Nmaxv(1)=Nmaxdy0;

%Getting out NmaxK and NmaxL values
for k = 2:Nmax+1
    NmaxK1 = Nmaxv(k-1);
    NmaxL1 = -(Nmaxy(k-1)^2-1)*Nmaxv(k-1)-Nmaxy(k-1);
    NmaxK2 = Nmaxv(k-1)+Nmaxh*NmaxL1;
    NmaxL2 = -((Nmaxy(k-1)+Nmaxh*NmaxK1)^2-1)*(Nmaxv(k-1)+Nmaxh*NmaxL1)-(Nmaxy(k-1)+Nmaxh*NmaxK1);
    NmaxK3 = Nmaxv(k-1)+Nmaxh*NmaxL1/4+Nmaxh*NmaxL2/4;
    NmaxL3 = -((Nmaxy(k-1)+Nmaxh*NmaxK1/4+Nmaxh*NmaxK2/4)^2-1)*(Nmaxv(k-1)+Nmaxh*NmaxL1/4+Nmaxh*NmaxL2/4)-(Nmaxy(k-1)+Nmaxh*NmaxK1/4+Nmaxh*NmaxK2/4);
  
    %Iterating over Nmaxy and Nmaxv using RK method
    Nmaxy(k)=Nmaxy(k-1)+(Nmaxh/6)*(NmaxK1+NmaxK2+NmaxK3);
    Nmaxv(k)=Nmaxv(k-1)+(Nmaxh/6)*(NmaxL1+NmaxL2+NmaxL3);
    Nmaxt=Nmaxt+Nmaxh;
    NmaxT=[NmaxT;Nmaxt];  

end
yNmax = Nmaxy(end); %error value we are comparing to


%Finding yN for the rest N=10,20,40,80,160

hvec=zeros(1,5);%our h values

for i=1:5
    %Initial values for each iteration
    t=0; tend=1; T=t;
    N=Nval(i);
    h=tend/N;
    hvec(i)=h;
    y0=1; y(1)=y0;
    dy0=0; v(1)=dy0;
    

for k=2:N+1
  
  %solving for K's and L's to use for RK method
  K1 = v(k-1);
  L1 = -(y(k-1)^2-1)*v(k-1)-y(k-1);
  K2 = v(k-1)+h*L1;
  L2 = -((y(k-1)+h*K1)^2-1)*(v(k-1)+h*L1)-(y(k-1)+h*K1);
  K3 = v(k-1)+h*L1/4+h*L2/4;
  L3 = -((y(k-1)+h*K1/4+h*K2/4)^2-1)*(v(k-1)+h*L1/4+h*L2/4)-(y(k-1)+h*K1/4+h*K2/4);
  
  %Iterating over y and v using RK method
  y(k)=y(k-1)+(h/6)*(K1+K2+K3);
  v(k)=v(k-1)+(h/6)*(L1+L2+L3);
    t=t+h;
    T=[T;t];
end
e(i)=abs(yNmax-y(end)); %getting the error
end
loglog(hvec, e,'o');
title('Loglog plot of error')
xlabel('h') 
ylabel('e')


%Note: order of accuracy looks like 2nd order from the graph based on how the order is decreasing



