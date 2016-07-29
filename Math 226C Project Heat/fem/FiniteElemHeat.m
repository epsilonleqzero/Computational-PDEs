%% Finite Element Method - Heat Equation
%
%   Written by Ted Kwan for Math 226C Project Heat
%
%   Returns u(x,t) to approximate the solution
%   to the Heat equation u_{t}-\Delta u=f
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
%
function [u,node,elem,t,Md,A,At]= FiniteElemHeat(node,elem,f,u0,t0,tf,dt,varargin)
    %% Initial Setup
    %
    t=[t0:dt:tf]'; Nt=length(t);
    N=length(node(:,1));
    [bdNode,~,isBdNode]=findboundary(elem);
    freeNode= find(~isBdNode);
    [A,M,area]=assemblematrix(node,elem);
    u=u0(node(:,1),node(:,2)); g_D=varargin{1};
    Md=spdiags(M,0,N,N); At=(Md+(dt/2).*A); Atr=((Md)-(dt/2).*A);
    for i=2:Nt
        uold=u; u=zeros(N,1);
        %% Boundary Condition.
        %
        u(bdNode) = g_D(node(bdNode,1),node(bdNode,2),t(i));
        %% Calculate RHS
        %
        ft=@(x,y) (f(x,y,t(i))+f(x,y,t(i-1)))/2;
        br=FiniteElemRHS(node,elem,ft,area,N);
        br=(dt).*br;  b1=Atr*uold;
        b=b1+br; r=b-At*u;
        %% Solve Au=r
        %
        u(freeNode)=At(freeNode,freeNode)\r(freeNode);
        %u=mg(At,r,elem);
    end
end