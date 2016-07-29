%% Adaptive Finite Element Method - Heat Equation
%
%   Written by Ted Kwan for Math 226C
%
%   Returns u(x,t) to approximate the solution
%   to the Heat equation u_{t}-\Delta u=f Using an
%   adaptive Crank-Nicolson scheme in both space and
%   time.
%
%%% Inputs
%
% * node - Nx2 Matrix of node coordinates.
% * elem - NTx3 Matrix of elements.
% * f - Function handle for the right hand side.
% * u0 - Function handle for the initial condition.
% * t0 - Initial time.
% * tf - Final time.
% * dt - Time step size.
% * varargin options:
% 1 function handle for g(x,t).
% 2 Type of adaptive refinement to use.
%
function [u,node,elem,t,M,A,At,errs]= astfemheat(node,elem,f,u0,t0,tf,dt,varargin)
    %% Initial Setup
    %
    maxitr=1000; t=zeros(maxitr,1);
    i=2; tcurr=t0; t(1)=tcurr;
    u=u0(node(:,1),node(:,2));
    g_D=varargin{1}; errs=cell(20,4); j=1;
    %% Main Loop
    %
    while(tcurr<tf && i<maxitr)
        uold=u; i=i+1;
        pde.f=@(x,y,t) f(x,y,t); pde.g_D=@ (x,y,z) g_D(x,y,z);
        %% Adaptive Solve
        %
        if(i>2)
            [u,M,A,At,node,elem,t(i),err,h1err,ns]= astfemsub1(node,elem,pde,uold,dt,t(i-1),tf,0);
        else
            [u,M,A,At,node,elem,t(i),err,h1err,ns]= astfemsub1(node,elem,pde,uold,dt,t(i-1),tf,1,u0);
        end
        %% Record Errors
        %
        errs{j,1}=err; errs{j,2}=ns; errs{j,3}=t(i);
        errs{j,4}=h1err; j=j+1;
    end
    %% Fixup for time
    %
    tswap=t(find(t)); t=[0;tswap];
end