%----------------------------------------------------------------
%C.2.1.3 Parameter Study of Solutions of an ODE System

%C.2.1.3.1 Problem 1: Particle Flow Past a Cylinder

%Going to use ode23 to solve this problem
R=2; 
yvec=[0.2,0.6,1.0,1.6];%contains our intital startig pts
for i=1:4
eqs = @(t,z) [1-((R.^2)*((z(1)).^2-(z(2)).^2))/((z(1)).^2+(z(2)).^2).^2;(-2*z(1)*z(2)*R.^2)/((z(1)).^2+(z(2)).^2).^2];%writing out eqns for ode23
[t,z] = ode23(@(t,z) eqs(t,z), [0,10], [-4;yvec(i)]);%ode23 being used to solve our ode
hold on
plot(z(:,1),z(:,2),'o');%grabbing all the points
title('Particle flow past a cylinder')
xlabel('x') 
ylabel('y')
axis equal;
end

%C.2.1.3.2 Problem 2: Motion of a particle


%Going to use ode23 to solve this problem

kvec=[0.020,0.065]; %different values for k
i=1; %counter
for k=kvec %We want to run twice(once for each k)
avec = [pi/6,pi/4,pi/3];%different angles
for a=avec
eqs = @(t,w) [w(3);w(4);-k*w(3)*sqrt((w(3).^2)+(w(4).^2));-9.81-k*abs(w(4))*sqrt((w(3).^2)+(w(4).^2))]; %setting up eqns for ode23
[t,w] = ode23(@(t,w) eqs(t,w), [0,2], [0;1.5;20*cos(a);20*sin(a)]); %I chose 2 for my t endpoint
figure(i)%new figure for each k
hold on
plot(w(:,1),w(:,2),'o');%want to grab all points
title('Motion of a particle')
xlabel('x') 
ylabel('y')
end
i=i+1;%increment
end