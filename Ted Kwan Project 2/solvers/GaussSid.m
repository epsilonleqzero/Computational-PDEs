%% Gauss-Sidel
%
%   Written by Ted Kwan for Math 226B
%   
%   This function implements the Gauss-Seidel
%   method using the code in the FDMcode 
%   document from the class notes by
%   Professor Chen.
function [u] = GaussSid(u,f)
    [n,m] = size(u);
    for j = 2:n-1
        for i = 2:m-1
            u(i,j)=(f(i,j)+u(i-1,j)+u(i+1,j)...
            +u(i,j-1)+u(i,j+1))/4;
        end
    end
end