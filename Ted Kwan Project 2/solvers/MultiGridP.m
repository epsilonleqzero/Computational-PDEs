%% Recursive Multigrid V-Cycle
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the recursive
%   multi-grid method for the Poisson
%   equation on the unit square with a
%   uniform grid.
function [u,hc] = MultiGridP(u,r,h,m,itr)
    %%% Coarse Grid
    %   Solve directly.
    if(itr==1)
        hc=h;
        u=Psolve(h,r);
        return ;
    end
    %%% Setup index map
    %
    Nc=(1/(2*h))+1;
    %%% Pre-Smoothing
    %
    for i=1:m
        u=GaussSid(u,r);
    end
    %%% Restriction
    %
    rf=(r-AuP(u)); rc=Res(rf,Nc);
    %%% Coarse Grid Correction
    %
    zc=zeros(Nc,Nc);
    [ec,hc]=MultiGridP(zc,rc,2*h,m,itr-1);
    %%% Post Smoothing
    %
    u=u+Prol(ec,Nc);
    for i=1:m
        u=GaussSidB(u,r);
    end
end