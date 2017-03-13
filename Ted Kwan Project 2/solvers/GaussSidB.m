%% Gauss-Sidel
%
%   Written by Ted Kwan for Math 226B
%   
%   This function implements the backwards
%   Gauss-Seidel method using the code in the
%   FDMcode document from the class notes
%   by Professor Chen Long.
function [u] = GaussSidB(u,f)
    [n,m] = size(u);
    for j=n-1:-1:2
        for i=m-1:-1:2
            u(i,j)=(f(i,j)+u(i-1,j)+u(i+1,j)+...
                   u(i,j-1)+u(i,j+1))/4;
        end
    end
end