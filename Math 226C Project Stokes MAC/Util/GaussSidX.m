%% Gauss-Sidel For Momentum Equation in U
%
%   Written by Ted Kwan for Math 226C Project 2.
%
%   This function does one iteration of gauss-seidel
%   to find the velocity function for u created
%   with the Marker and Cell (MAC) scheme.

function [u] = GaussSidX(u,p,f,omega)
[n,m] = size(p);
%% Boundary
%
j = 2:n;
u(1,j)=(1-omega)*u(1,j)+omega*(f(1,j)+(4/3)*u(2,j)+u(1,j-1)+...
    u(1,j+1)-(p(1,j)-p(1,j-1)))/6; 
u(n,j)=(1-omega)*u(n,j)+omega*(f(n,j)+(4/3)*u(n-1,j)+u(n,j-1)+...
    u(n,j+1)-(p(n,j)-p(n,j-1)))/6;
%% Interior
%
% case 1 (red points): mod(i+j,2) == 0
i = 2:2:m-1; j = 2:2:n;
u(i,j) = (1-omega)*u(i,j)+omega*(f(i,j) + u(i-1,j) + u(i+1,j) ... 
        + u(i,j-1) + u(i,j+1)-(p(i,j)-p(i,j-1)))/4;
i = 3:2:m-1; j = 3:2:n;
u(i,j) = (1-omega)*u(i,j)+omega*(f(i,j) + u(i-1,j) + u(i+1,j) ... 
        + u(i,j-1) + u(i,j+1)-(p(i,j)-p(i,j-1)))/4;
% case 2 (black points): mod(i+j,2) == 1
i = 2:2:m-1; j = 3:2:n;
u(i,j) = (1-omega)*u(i,j)+omega*(f(i,j) + u(i-1,j) + u(i+1,j) ... 
        + u(i,j-1) + u(i,j+1)-(p(i,j)-p(i,j-1)))/4;
i = 3:2:m-1; j = 2:2:n;
u(i,j) = (1-omega)*u(i,j)+omega*(f(i,j) + u(i-1,j) + u(i+1,j) ... 
        + u(i,j-1) + u(i,j+1)-(p(i,j)-p(i,j-1)))/4;
end