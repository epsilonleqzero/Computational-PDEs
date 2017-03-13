%% Gauss-Sidel Red Black
%
%   Obtained from FDMcode and used for
%   Project 2
function [u] = GaussSidRB(u,f)
[n,m] = size(u);
% case 1 (red points): mod(i+j,2) == 0
i = 2:2:m-1; j = 2:2:n-1;
u(i,j) = (f(i,j) + u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
i = 3:2:m-1; j = 3:2:n-1;
u(i,j) = (f(i,j) + u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
% case 2 (black points): mod(i+j,2) == 1
i = 2:2:m-1; j = 3:2:n-1;
u(i,j) = (f(i,j) + u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
i = 3:2:m-1; j = 2:2:n-1;
u(i,j) = (f(i,j) + u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;

end