%% Gauss-Sidel For Momentum Equation in V
%
%   Written by Ted Kwan for Math 226C Project 2.
%
%   This function does one iteration of gauss-seidel
%   to find the velocity function for v created
%   with the Marker and Cell (MAC) scheme.

function [v] = GaussSidY(v,p,f,omega)
[n,m] = size(p);
%% Boundary
%
i = 2:m;
v(i,1)=(1-omega)*v(i,1)+omega*(f(i,1)+(4/3)*v(i,2)+v(i-1,1)+...
    v(i+1,1)-(p(i-1,1)-p(i,1)))/6;
v(i,m)=(1-omega)*v(i,m)+omega*(f(i,m)+(4/3)*v(i,m-1)+v(i-1,m)+...
    v(i+1,m)-(p(i-1,m)-p(i,m)))/6;
%% Interior
%
% case 1 (red points): mod(i+j,2) == 0
i = 2:2:m; j = 2:2:n-1;
v(i,j) = (1-omega)*v(i,j)+omega*(f(i,j) + v(i-1,j) + v(i+1,j) ... 
        + v(i,j-1) + v(i,j+1)-(p(i-1,j)-p(i,j)))/4;
i = 3:2:m; j = 3:2:n-1;
v(i,j) = (1-omega)*v(i,j)+omega*(f(i,j) + v(i-1,j) + v(i+1,j) ... 
        + v(i,j-1) + v(i,j+1)-(p(i-1,j)-p(i,j)))/4;
% case 2 (black points): mod(i+j,2) == 1
i = 2:2:m; j = 3:2:n-1;
v(i,j) = (1-omega)*v(i,j)+omega*(f(i,j) + v(i-1,j) + v(i+1,j) ... 
        + v(i,j-1) + v(i,j+1)-(p(i-1,j)-p(i,j)))/4;
i = 3:2:m; j = 2:2:n-1;
v(i,j) = (1-omega)*v(i,j)+omega*(f(i,j) + v(i-1,j) + v(i+1,j) ... 
        + v(i,j-1) + v(i,j+1)-(p(i-1,j)-p(i,j)))/4;
end