%% Distributed Gauss-Seidel for Maxwell Equations
%
%   Written by Ted Kwan for Math 226C Project Maxwell
%
%   This function runs the distributed Gauss-Seidel
%   method residual correction method.
%
function [eu,ep] = DGSMaxwell(ru,rp,ep,B,R,Rp)
    %% One Step of G-S on u
    %
    eu=R\(ru);
    %% Poisson equation for dq
    %
    b=-B*eu+rp;
    b=b-mean(b);
    dq=Rp\(b);
    %% Update Corrections
    %
    eu=eu+B'*dq;
    ep=ep-(B*B')*dq;
end