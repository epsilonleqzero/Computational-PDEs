%% Two-Grid Method
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the two-grid
%   method for the Poisson equation on the
%   unit square with a uniform grid.
function [u] = TwoGridP(u,f,h,m,itr)
    %%% Setup Index maps
    %
    Nc=(1/(2*h))+1;
    %%% Pre-Smoothing
    %
    for i=1:m
        u=GaussSid(u,f);
    end
    %%% Solve at Coarse-Grid
    %
    Au=AuP(u); rc=Res(f-Au,Nc);
    ec=Psolve(2*h,rc);
    %%% Post Smoothing
    %
    u=u+Prol(ec,Nc);
    for i=1:m
        u=GaussSidB(u,f);
    end
end