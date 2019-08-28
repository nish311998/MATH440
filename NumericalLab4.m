% Computational Lab 4
% By: Nishanth Gona
% May 3, 2019


%Part 1.1
h=0.2;N=19;M=11;%made M = 11 to have ghost points, and to have h = 0.2
psize = N*M + 2*(N) + 2*M + 4; %getting our problem size
A = sparse(psize, psize);
b=zeros(psize, 1);


rowind = 1; 

delx = 0.2;%initialize delx to 0.2
dely = 0.2;%initialize dely to 0.2

%building the main part of my matrix 
for j = 2:M+1
    for i = 2:N+1
        cind1 = index(i,j,N);
        cind2 = index(i,j+1,N);
        cind3 = index(i,j-1,N);
        cind4 = index(i+1,j,N);
        cind5 = index(i-1,j,N);
        A(rowind, cind1) = -2*(1+(delx/dely).^2);
        A(rowind, cind2) = (delx/dely).^2;
        A(rowind, cind3) = (delx/dely).^2;
        A(rowind, cind4) = 1;
        A(rowind, cind5) = 1;
        b(rowind) = 0;
        rowind  = rowind + 1;
    end
end

%building bc 1
%rowind = 1;  
for j = 2:M+1
    i=1;
    colind = index(i,j,N);
    A(rowind, colind) = 1;
    b(rowind) = 300;
    rowind  = rowind + 1;
end

%building bc 2
%rowind = 1;
for j = 2:M+1
    i=N+2;
    colind = index(i,j,N);
    A(rowind, colind) = 1;
    b(rowind) = 600;
    rowind  = rowind + 1;
end

%building bc 3
%rowind = 1;
for i = 2:N+1 
    jtemp = 2;
    colind1 = index(i,jtemp+1,N);
    colind2 = index(i,jtemp-1,N);
    A(rowind, colind1) = 1;
    A(rowind, colind2) = -1;
    b(rowind) = 0;
    rowind  = rowind + 1;
end

%building bc 4
%rowind = 1;
for i = 2:N+1 
    jtemp = M+1;
    colind1 = index(i,jtemp+1,N);
    colind2 = index(i,jtemp-1,N);
    A(rowind, colind1) = 1;
    A(rowind, colind2) = -1;
    b(rowind) = 0;
    rowind  = rowind + 1;
end


%Corner pts
colind1 = index(1,1,N);
A(rowind,colind1) = 1;
b(rowind)=0;
rowind = rowind+1;

colind2 = index(N+2,1,N);
A(rowind,colind2) = 1;
b(rowind) = 0;
rowind = rowind+1;

colind3 = index(1,M+2,N);
A(rowind,colind3) = 1;
b(rowind) = 0;
rowind = rowind+1;

colind4 = index(N+2,M+2,N);
A(rowind,colind4) = 1;
b(rowind)=0;



%Solution
y = A\b;    


%Turning vector into matrix
resultmatrix = zeros(M, N+2);
for j = 2:M+1
    for i = 1:N+2
        resultmatrix(j-1,i) = y(index(i,j,N));
    end
end

xint = linspace(0, 4, 21);
yint = linspace(0, 2, 11);
[X, Y] = meshgrid(xint,yint);

figure(1)
surf(X, Y, resultmatrix);
xlabel('x');
ylabel('y');
zlabel('Temperature');
title('Temperature Distribution');

%To display using imagesc
figure(2)
imagesc(0:0.2:4,0:0.2:2,resultmatrix);
xlabel('x');
ylabel('y');
title('Temperature Distribution in color');



%Values at T(2,1)
pt = index(11, 7, N);
ypart1 = y(pt);


%Proving the true solution

%After finding the soln T(x,y) = 75x + 300

%analytic matrix
provematrix = zeros(M, N+2);
for j = 2:M+1
    for i = 1:N+2
        provematrix(j-1,i) = 75*xint(i) + 300;
    end
end

%Getting error
finalmatrix = abs(provematrix - resultmatrix);

%error
figure(3)
mesh(X, Y, finalmatrix);
xlabel('x');
ylabel('y');
zlabel('Error');
title('Error of Temperature Distribution between Analytic Soln and Numerical Soln');

%Now using LU Factorization to solve the same problem, part 2 of lab

[L,U] = lu(A);

y2 = inv(L).*b;  
x = inv(U).*y2;
%%

%Part 1.2
h=0.1;N=39;M=21;%made M = 21 to have ghost points, and to have h = 0.2
psize = N*M + 2*N + 2*M +4;
A = sparse(psize, psize);
b=zeros(psize, 1);


rowind = 1; 

delx = 0.1;%initialize delx to 0.1
dely = 0.1;%initialize dely to 0.1

%building the main part of my matrix 
for j = 2:M+1
    for i = 2:N+1
        cind1 = index(i,j,N);
        cind2 = index(i,j+1,N);
        cind3 = index(i,j-1,N);
        cind4 = index(i+1,j,N);
        cind5 = index(i-1,j,N);
        A(rowind, cind1) = -2*(1+(delx/dely).^2);
        A(rowind, cind2) = (delx/dely).^2;
        A(rowind, cind3) = (delx/dely).^2;
        A(rowind, cind4) = 1;
        A(rowind, cind5) = 1;
        b(rowind) = 0;
        rowind  = rowind + 1;

    end
end

%building bc 1
%rowind = 1;
for j = 2:M+1
    i=1;
    colind = index(i,j,N);
    A(rowind, colind) = 1;
    b(rowind) = 300;
    rowind  = rowind + 1;
end

%building bc 2
%rowind = 1;
for j = 2:M+1
    i=N+2;
    colind = index(i,j,N);
    A(rowind, colind) = 1;
    b(rowind) = 600;
    rowind  = rowind + 1;
end

%building bc 3
%rowind = 1;
for i = 2:N+1
    colind1 = index(i,3,N);
    colind2 = index(i,1,N);
    A(rowind, colind1) = 1;
    A(rowind, colind2) = -1;
    b(rowind) = 0;
    rowind  = rowind + 1;
end

%building bc 4
%rowind = 1;
for i = 2:N+1
    colind1 = index(i,M+2,N);
    colind2 = index(i,M,N);
    A(rowind, colind1) = 1;
    A(rowind, colind2) = -1;
    b(rowind) = 0;
    rowind  = rowind + 1;
end

%Corner pts
colind1 = index(1,1,N);
A(rowind,colind1) = 1;
b(rowind)=0;
rowind = rowind+1;

colind2 = index(N+2,1,N);
A(rowind,colind2) = 1;
b(rowind) = 0;
rowind = rowind+1;

colind3 = index(1,M+2,N);
A(rowind,colind3) = 1;
b(rowind) = 0;
rowind = rowind+1;

colind4 = index(N+2,M+2,N);
A(rowind,colind4) = 1;
b(rowind)=0;

%Solution
y = A\b;    

%Turning vector back to matrix
resultmatrix2 = zeros(M, N+2);
for j = 2:M+1
    for i = 1:N+2
        resultmatrix2(j-1,i) = y(index(i,j,N));
    end
end

xint = linspace(0, 4, 41);
yint = linspace(0, 2, 21);
[X, Y] = meshgrid(xint,yint);
figure(4)
surf(X, Y, resultmatrix2);
xlabel('x');
ylabel('y');
zlabel('Temperature');
title('Temperature Distribution');

%To display using imagesc
figure(5)
imagesc(0:0.2:4,0:0.2:2, resultmatrix2);
xlabel('x');
ylabel('y');
title('Temperature Distribution in color');



%Values at T(2,1)
pt = index(21, 12, N);
ypart2 = y(pt);


%analytic matrix
provematrix2 = zeros(M, N+2);
for j = 2:M+1
    for i = 1:N+2
        provematrix2(j-1,i) = 75*xint(i) + 300;
    end
end

%Getting error
finalmatrix2 = abs(provematrix2 - resultmatrix2);

%error
figure(6)
mesh(X, Y, finalmatrix2);
xlabel('x');
ylabel('y');
zlabel('Error');
title('Error of Temperature Distribution between Analytic Soln and Numerical Soln');

%Now using LU Factorization to solve the same problem

[L,U] = lu(A);

y2 = inv(L).*b;
x = inv(U).*y2;

fprintf('T(2,1) for the first part is %i\n', + ypart1);
fprintf('T(2,1) for the second part is %i\n', + ypart2);
fprintf('Ansatz for the value at T(2,1) would be 450\n');
fprintf('This is proven by our error plots in Figure 3 and 6, as seen by the very minimal error\n');


%helper function for index
function val = index(i, j, N)
val = (j-1).*(N+2)+i;
end


