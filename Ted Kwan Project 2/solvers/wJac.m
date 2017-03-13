%% Weighted Jacobi
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the weighted
%   Jacobi method.
function [u] = wJac(u,f,omega)
    [n,m]=size(u);i=2:m-1;j=2:n-1; uJ=u;
    %%% Normal iteration
    %
    uJ(i,j) =(f(i,j)+u(i-1,j)+u(i+1,j) ...
            +u(i,j-1)+u(i,j+1))/4;
    %%% Weighted correction
    %
    u=omega*u+(1-omega)*uJ;
end