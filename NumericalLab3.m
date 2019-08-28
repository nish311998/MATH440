%--------------------------------------------------------------------------
% Computational Lab 3
% By: Nishanth Gona
% April 9, 2019

%Note:Please run separetly for each part when you see '%%' this means a
%separate section
%% 

%Part 3

%First make our big U vector

%Part 3

%Building all at once

%stable

a=1;
tn=10000;%time steps
delt = 2/tn; %tau=2
n = 20; %lets choose 20
L=1; %end value
h = L/(n+1);  %Grid spacing 
u = zeros(tn,n+2);%want matrix of this size
u(1:tn/2,1) = 1;%until t=1 vals=1
u(:,n+2) = 0;
k=a*delt/(h.^2)%at/h^2
    for j = 1:tn
        for i = 2:n+1   
            u(j+1, i) =k*u(j,i+1)+(1-2*k)*u(j,i)+k*(u(j,i-1));%center difference
        end
        u(j+1,n+2) = 2*k*u(j,n+1)+(1-2*k)*u(j,n+2);%boundary condition
    end

%This is the surface plot
[H,T] = meshgrid(0:h:1,0:delt:2);

mesh(H,T,u);
xlabel('space');
ylabel('time');
zlabel('u');
title('Temp in rod Stable');

fprintf('h is %i\n', + h);
fprintf('t is %i\n', + delt);
fprintf('k is %i\n', + k);

%% 

%unstable

a=1;
tn=1500;%time steps
delt = 2/tn; %tau=2
n = 20; %lets choose 20
L=1; %end value
h = L/(n+1);  %Grid spacing 
u = zeros(tn,n+2);%want matrix of this size
u(1:tn/2,1) = 1;%until t=1 vals=1
u(:,n+2) = 0;
k=a*delt/(h.^2);%at/h^2
    for j = 1:tn
        for i = 2:n+1   
            u(j+1, i) =k*u(j,i+1)+(1-2*k)*u(j,i)+k*(u(j,i-1));%center difference
        end
        u(j+1,n+2) = 2*k*u(j,n+1)+(1-2*k)*u(j,n+2);%boundary condition
    end

%This is the surface plot
[H,T] = meshgrid(0:h:1,0:delt:2);

mesh(H,T,u);
xlabel('space');
ylabel('time');
zlabel('u');
title('Temp in rod Unstable');

fprintf('h is %i\n', + h);
fprintf('t is %i\n', + delt);
fprintf('k is %i\n', + k);
%% 

%Part 4
%ode23 and %ode23s

%ode23
a=1;
n=10;
L=1;
tn=10000;
delt = 2/tn;
data=zeros(3);

for i = 1:3 %looping 3 times
    start = zeros(n+2,1);

    h=L/n+1;  
    k=a*delt/(h.^2);
    rel=10.^-6; %these are our relative tolerances
    options = odeset('RelTol',rel); %setting our RelTol
    tic;
    [t,u] = ode23(@(t,u)ode(t,u,n,k,h), [0,2], start , options); %ode23 calling ode function
    time=toc;

    %making chart of values
    data(i,1) = length(t);
    data(i,2) = time;
    data(i,3) = max(diff(t));
    n=n*2;
end

fprintf('Time steps for n=10 is %i\n', data(1,1)-1);
fprintf('Time steps for n=20 is %i\n', data(2,1)-1);
fprintf('Time steps for n=40 is %i\n', data(3,1)-1);

fprintf('Elapsed CPU time for n=10 is %i secs\n', data(1,2));
fprintf('Elapsed CPU time for n=20 is %i secs\n', data(2,2));
fprintf('Elapsed CPU time for n=40 is %i secs\n', data(3,2));

fprintf('Max time step for n=10 is %i\n', data(1,3));
fprintf('Max time step for n=20 is %i\n', data(2,3));
fprintf('Max tiem step for n=40 is %i\n', data(3,3));

%check data matrix for values 
%% 


%ode23s
a=1;
n=10;
L=1;
tn=10000;
delt = 2/tn;
data=zeros(3);

for i = 1:3
    
    start = zeros(n+2,1);
    h=L/n+1;  
    k=a*delt/(h.^2);
    rel=10.^-6; %these are our relative tolerances
    options = odeset('RelTol',rel); %setting our RelTol
    tic;
    [t,u] = ode23s(@(t,u)ode(t,u,n,k,h), [0,2], start , options); %ode23 calling the ode function
    time=toc;
    
    %Chart values
    data(i,1) = length(t);
    data(i,2) = time;
    data(i,3) = max(diff(t));
    n=n*2;
end
fprintf('Time steps for n=10 is %i\n', data(1,1)-1);
fprintf('Time steps for n=20 is %i\n', data(2,1)-1);
fprintf('Time steps for n=40 is %i\n', data(3,1)-1);

fprintf('Elapsed CPU time for n=10 is %i secs\n', data(1,2));
fprintf('Elapsed CPU time for n=20 is %i secs\n', data(2,2));
fprintf('Elapsed CPU time for n=40 is %i secs\n', data(3,2));

fprintf('Max time step for n=10 is %i\n', data(1,3));
fprintf('Max time step for n=20 is %i\n', data(2,3));
fprintf('Max tiem step for n=40 is %i\n', data(3,3));

%% 


%Part 5
%Use jacobian and jpattern

%a
%ode23
a=1;
n=10;
L=1;
tn=11;
delt = 2/tn;
start = zeros(n+2,1);

    h=L/n+1;  
    k=a*delt/(h.^2);
    u = zeros(tn,n+2);
    u(1:tn/2,1) = 1;
    u(:,n+2) = 0;
    %want to create u before to use with Jacobian
    for j = 1:tn
        for i = 2:n+1   
            u(j+1, i) =k*u(j,i+1)+(1-2*k)*u(j,i)+k*(u(j,i-1));
        end
        u(j+1,n+2) = 2*k*u(j,n+1)+(1-2*k)*u(j,n+2);
    end
    rel=10.^-6; %these are our relative tolerances
    options = odeset('Jacobian',u,'RelTol',rel); 
    %options = odeset('RelTol',rel); %setting our RelTol
    tic;
    help = @(t,u)ode(t,u,n,k,h);
    [t,u] = ode23s(help, [0,2], start , options);  
    time = toc;
fprintf('Elapsed CPU time for n=10 is %i secs\n', time);

%% 

%b
a=1;
n=10;
L=1;
tn=10000;
delt = 2/tn;
start = zeros(n+2,1);

    h=L/n+1;  
    k=a*delt/(h.^2);
    u = zeros(tn,n+2);
    u(1:tn/2,1) = 1;
    u(:,n+2) = 0;
    %want to create u before to use with Jpattern
    for j = 1:tn
        for i = 2:n+1   
            u(j+1, i) =k*u(j,i+1)+(1-2*k)*u(j,i)+k*(u(j,i-1));
        end
        u(j+1,n+2) = 2*k*u(j,n+1)+(1-2*k)*u(j,n+2);
    end
    rel=10.^-6; %these are our relative tolerances
    %we want sparse
    %options = odeset('Jpattern',u,'RelTol',rel); %wanted to use this but
    %got an array out of bounds error
    options = odeset('RelTol',rel); 
    tic;
    [t,u] = ode23s(@(t,u)ode(t,u,n,k,h), [0,2], start, options); 
    time = toc;
fprintf('Elapsed CPU time for n=10 is %i secs\n', time);


%Also no search field named odefile to do the third option
%% 
%Part 6
a=1;
tn=10000;
delt = 2/tn; %changing
n = 20; %lets choose 20
L=1; %end value
h = L/(n+1);  %Grid spacing 
u = zeros(tn,n+2);
u(1:tn/2,1) = 1;
u(:,n+2) = 0;
k=a*delt/(h.^2)
%this is making u again
    for j = 1:tn
        for i = 2:n+1   
            u(j+1, i) =k*u(j,i+1)+(1-2*k)*u(j,i)+k*(u(j,i-1));
        end
        u(j+1,n+2) = 2*k*u(j,n+1)+(1-2*k)*u(j,n+2);
    end
   %obtaining 4 row vectors
 uvec=[u(2500,:); u(5000,:); u(7500,:); u(10000,:)];
 x=linspace(0,1,n+2);
 figure(1)
 plot(x,uvec)
 xlabel('x');
 ylabel('u')
 title('Analysis at t=0.5,1,1.5,2');
 [H,T] = meshgrid(0:h:1,0:delt:2);

 
 %Pretty much same graph we had in part 1
figure(2);
mesh(H,T,u);
xlabel('space');
ylabel('time');
zlabel('u');
title('Temp in rod');


%Answers to problem 6:


%1.Use a stiff method because it will reduce the number of time steps

%2. The jacobian matrix is quite large so it cannot be a band matrix, and
%takes the shape of a spare matrix

%3.For the implicit ode, the solver for sparse matrix should be used as it
%reduces the number of calculations we need to perform
%% 
%function we are using for ode solvers
function out = ode(t, u, n, k, h)
out = zeros(n+2,1);
out(1)=(t<1);
for i = 2:n+1   
     out(i)=k*u(i+1)+(1-2*k)*u(i)+k*(u(i-1));
end
out(n+2) = 2*k*u(n+1)+(1-2*k)*u(n+2);

end
