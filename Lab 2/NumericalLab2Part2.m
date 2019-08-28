%---------------------------------------------------------------
%C2.1.2 Stability Investigation of a RK Method

%C2.1.2.1 Constant Stepsize Experiment
Nval=[125, 250, 500, 1000, 2000]; %different N
k1=0.04;k2=10.^4;k3=3*10.^7; %different k

for i = 1:5
    %Initial Values for each iteration
    t=0; tend=1; T=t;
    N=Nval(i);
    h=tend/N;
    x1(1)=1; x2(1)=0; x3(1)=0;
for k=2:N+1
  
  %solving for K's, L's, and J's iteratively using Rk from part 1
  K1 = -k1*x1(k-1)+k2*(x2(k-1))*(x3(k-1));
  L1 = k1*x1(k-1)-k2*((x2(k-1))*(x3(k-1)))-k3*((x2(k-1)).^2);
  J1 = k3*((x2(k-1)).^2);
  K2 = -k1*(x1(k-1)+h*K1)+k2*(x2(k-1)+h*L1)*(x3(k-1)+h*J1);
  L2 = k1*(x1(k-1)+h*K1)-k2*((x2(k-1)+h*L1)*(x3(k-1)+h*J1))-k3*((x2(k-1)+h*L1).^2);
  J2 = k3*((x2(k-1)+h*L1).^2);
  K3 = -k1*(x1(k-1)+h*K1/4+h*K2/4)+k2*(x2(k-1)+h*L1/4+h*L2/4)*(x3(k-1)+h*J1/4+h*J2/4);
  L3 = k1*(x1(k-1)+h*K1/4+h*K2/4)-k2*((x2(k-1)+h*L1/4+h*L2/4)*(x3(k-1)+h*J1/4+h*J2/4))-k3*((x2(k-1)+h*L1/4+h*L2/4).^2);
  J3 = k3*((x2(k-1)+h*L1/4+h*L2/4).^2);
  
  %Iterating over x1,x2, x3 and updating values
  x1(k)=x1(k-1)+(h/6)*(K1+K2+K3);
  x2(k)=x2(k-1)+(h/6)*(L1+L2+L3);
  x3(k)=x3(k-1)+(h/6)*(J1+J2+J3);
  t=t+h;
  T=[T;t];
    
end

hold on
figure(i);
loglog(T,x1,'o');
loglog(T,x2,'o');
loglog(T,x3,'o');
title('Robertsons Problem')
xlabel('T') 
ylabel('xvector')
end

%Note the smallest N at which the points converge is 1000, which is the 4th
%figure


%C.2.1.2.2 Adaptive Stepsize Experiment Using Matlab Functions
%part 1

relvec=[10.^-3,10.^-4,10.^-5,10.^-6]; %these are our relative tolerances
stepvec=zeros(4);
i=1;
for rel=relvec %looping through all elements in relveci=1;
options = odeset('RelTol',rel); %setting our RelTol

eqs = @(t,x) [-k1*x(1)+k2*x(2)*x(3);k1*x(1)-k2*x(2)*x(3)-k3*(x(2)).^2;k3*(x(2)).^2]; %setting up equations for ode23
[t,x] = ode23(@(t,x) eqs(t,x), [0,1], [1;0;0], options); %ode23
figure(6);
hold on
plot(t(2:end),diff(t));%want vectors to match up so we start from second element in t
title('Robertsons Problem Nonstiff')
xlabel('t') 
ylabel('h')
stepvec(i)=length(t);
i=i+1;
end

%As seen in stepvec our steps for each tolerance using ode23 are obtained by our t
%values: so we have 10^-3:867, 10^-4:868, 10^-5:869, 10^-6:869

%part2
stepvec2=zeros(4);
i=1;
for rel=relvec

options = odeset('RelTol',rel); %setting out RelTol

eqs = @(t,x) [-k1*x(1)+k2*x(2)*x(3);k1*x(1)-k2*x(2)*x(3)-k3*(x(2)).^2;k3*(x(2)).^2]; %setting up equations for ode23s
[t,x] = ode23s(@(t,x) eqs(t,x), [0,1000], [1;0;0], options);%ode23s
figure(7);
hold on
plot(t(2:end),diff(t)); %want vectors to match up so we start from second element in t
title('Robertsons Problem Stiff')
xlabel('t') 
ylabel('h')
stepvec2(i)=length(t);
i=i+1;
end

%As seen in stepvec2 our steps for each tolerance using ode23s are obtained by our t
%values: so we have 10^-3:31, 10^-4:38, 10^-5:49, 10^-6:62