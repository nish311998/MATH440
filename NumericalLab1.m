%--------------------------------------------------------------------------
% Computational Lab 1
% By: Nishanth Gona
% Feb. 17, 2019
%
% Solves the following problem:
% -kappa*u_ZZ + v*p*c(u_Z) = Q(z)
% u(0) = Intial T, u(L) = -kappa*T_Z(L)= k*(T(L)-Tout)
%
function lab1
    % constants in the lab
    L=10;%end of interval
    aa=1;%part of piecewise
    bb=3;%part of piecewise
    Qinit = 50;
    kappa = 0.5;
    k=10;
    p=1;
    c=1;
    Tout=300;
    Tinit=400;
    v=0;
    
    %First Chart
    n = 10; %our initial N 
    
    for i = 1:4 %using for n=10,20,40,80
    x = (0:(n))'*L/(n); %our x values in the grid
    h = L/(n);  %Grid spacing
    f = zeros(n,1);      %External forcing 
    %bc_0 = Tinit;   %Change these values to implement different Dirichlet
    %bc_1 = (6000*h+b(end-1))/(20*h +1);   %Nuemann boundary conditions  
    
    
    %Build matrix   A = second derivative matrix
    %               B = center difference
    A = zeros(n+1); %Exclude j = 0 and j = n+1
    B = zeros(n+1);
    for j = 2:n
        A(j, j) = -2;
        A(j, j+1) = 1;
        A(j, j-1) = 1;
        
        B(j, j+1) = 1/2*h;
        B(j, j-1) = -1/2*h;        
    end
    A(1, 1) = -1; %Initial Boundary Condition
    A(n+1, n+1) = -20*h - 1; A(n+1, n) = 1; %End Boundary condition
    B(1, 1) = 1; %Initial Boundary Condition
    B(n+1, n+1) = 20*h + 1; B(n+1, n) = -1;%End Boundary condition
    %Build the "right hand side"
    %We need n-2 intermediate points plus two boundary conditions here
    %building our f
      for z = 1:length(x)  
        i=x(z);
        if i <= aa
            f(z) = 0;
        elseif i >= aa && i<=bb
            f(z) = Qinit*sin((i-aa)*pi/(bb-aa));%forcing 
        else
            f(z)=0;
        end
      end
    
    %creating right hand side
    b = f*h^2;
    b(1) = Tinit;
    b(end) = (6000*h+b(end-1))/(20*h +1); 
    
    M = -kappa*A + v*p*c*B;
    M(1,1)=1;
    M(n+1,n+1) = 20*h+1; M(n+1,n) = -1;
    %Solve By = b using MATLABs "Backslash" operator.
    y = M\b;    %Solution. 
    figure(1)
    title('Varying N')
    xlabel('Z') 
    ylabel('T(z)')
    %legend('N=10', 'N=20', 'N=40', 'N=80')  Gave me an error when I tried
    %doing this

    %plot(x, y, 'o', 0, bc_0, 'o', L, bc_1, 'o', 'linewidth', 3)
    plot(x,y,'o')
    n = n*2; %Get our new n for the next iteration
    hold on
    end
    hold off
    
    %Second Chart
    diffv=[0.1;0.5;1;10]; %all of our diff v
    n=40;
    for i = 1:length(diffv)
    v = diffv(i);
    x = (0:(n))'*L/(n);
    h = L/(n);  %Grid spacing
    f = zeros(n,1);      %External forcing 
    
    %Qinit*sin((z-a)pi/(b-a));
    %bc_0 = Tinit;   %Change these values to implement different Dirichlet
    %bc_1 = (6000*h+b(end-1))/(20*h +1);   %Nuemann boundary conditions  IDKKKKKKK
    
    %Build matrix   A = second derivative matrix
    %               B = center difference
    A = zeros(n+1); %Exclude j = 0 and j = n+1
    B = zeros(n+1);
    for j = 2:n
        A(j, j) = -2;
        A(j, j+1) = 1;
        A(j, j-1) = 1;
        
        B(j, j+1) = 1/2*h;
        B(j, j-1) = -1/2*h;        
    end
    A(1, 1) = -1;
    A(n+1, n+1) = -20*h - 1; A(n+1, n) = 1;
    B(1, 1) = 1; 
    B(n+1, n+1) = 20*h + 1; B(n+1, n) = -1;
    %Build the "right hand side"
    %We need n-2 intermediate points plus two boundary conditions here
    %building our f
      for z = 1:length(x)  
        i=x(z);
        if i <= aa
            f(z) = 0;
        elseif i >= aa && i<=bb
            f(z) = Qinit*sin((i-aa)*pi/(bb-aa));
        else
            f(z)=0;
        end
      end
    
    %creating right hand side
    b = f*h^2;
    b(1) = Tinit;
    b(end) = (6000*h+b(end-1))/(20*h +1);     
    M = -kappa*A + v*p*c*B;
    M(1,1)=1;
    M(n+1,n+1) = 20*h+1; M(n+1,n) = -1;
    %Solve By = b using MATLABs "Backslash" operator.
    y = M\b;    %Solution. 
    figure(2)
    xlabel('Z') 
    ylabel('T(z)')
    title('Varying v')

    %plot(x, y, 'o', 0, bc_0, 'o', L, bc_1, 'o', 'linewidth', 3)
    plot(x,y,'o')
    hold on
    end
    hold off


    % v = 10 and different N
    diffn=[10;20;40];
    v=10;
    for i = 1:length(diffn)
    n = diffn(i);
    x = (0:(n))'*L/(n);
    h = L/(n);  %Grid spacing
    f = zeros(n,1);      %External forcing 

  
    
    %Qinit*sin((z-a)pi/(b-a));
    %bc_0 = Tinit;   %Change these values to implement different Dirichlet
    %bc_1 = (6000*h+b(end-1))/(20*h +1);   %Nuemann boundary conditions  
    
    %Build matrix   A = second derivative matrix
    %               B = center difference
    A = zeros(n+1); %Exclude j = 0 and j = n+1
    B = zeros(n+1);
    for j = 2:n
        A(j, j) = -2;
        A(j, j+1) = 1;
        A(j, j-1) = 1;
        
        B(j, j+1) = 1/2*h;
        B(j, j-1) = -1/2*h;        
    end
    A(1, 1) = -1;
    A(n+1, n+1) = -20*h - 1; A(n+1, n) = 1;
    B(1, 1) = 1; 
    B(n+1, n+1) = 20*h + 1; B(n+1, n) = -1;
    %Build the "right hand side"
    %We need n-2 intermediate points plus two boundary conditions here
    %building our f
      for z = 1:length(x)  
        i=x(z);
        if i <= aa
            f(z) = 0;
        elseif (i >= aa) & (i <= bb)
            f(z) = Qinit*sin((i-aa)*pi/(bb-aa));
        else
            f(z)=0;
        end
      end
    
    %Creating right hand side
    b = f*h^2;
    b(1) = Tinit;
    b(end) = (6000*h+b(end-1))/(20*h +1);     
    
    
    M = -kappa*A + v*p*c*B;
    M(1,1)=1;
    M(n+1,n+1) = 20*h+1; M(n+1,n) = -1;
    
    %Solve By = b using MATLABs "Backslash" operator.
    y = M\b;    %Solution. 
    title('Varying N and v=10')
    xlabel('Z') 
    ylabel('T(z)')
    figure(3)
    %plot(x, y, 'o', 0, bc_0, 'o', L, bc_1, 'o', 'linewidth', 3)
    plot(x,y,'o')
    hold on
    end
    hold off
    
    
end


%Regarding how to get rid of spurious oscillations:
%The best way to get rid of spurious oscillations is to have a smaller h
%and we can do this by increasing the number of points
