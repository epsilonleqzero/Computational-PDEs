%% AFEM Heat Solving method
%
%   Written by Ted Kwan for Math 226C Project Heat.
%
%   This function sets up and solves the heat equation on
%   a single time step.
%
function [u,Md,A,At] = afemheatsolve(node,elem,g_D,f,uold,dt)
    %% Initial Setup
    %
    [A,M,area]=assemblematrix(node,elem);
    [bdNode,~,isBdNode]=findboundary(elem);
    freeNode= find(~isBdNode);
    N=length(node(:,1)); Md=spdiags(M,0,N,N);
    At=(Md+(dt/2).*A); Atr=((Md)-(dt/2).*A);
    u=zeros(N,1);
    %% Dirichlet Boundary Condition
    %
    u(bdNode) = g_D(node(bdNode,1),node(bdNode,2));
    br=FiniteElemRHS(node,elem,f,area,N);
    %% Create RHS
    %
    br=(dt).*br; b1=Atr*uold;
    b=b1+br; r=b-At*u;
    %%% End RHS
    %% Solve Au=r
    %
    u(freeNode)=At(freeNode,freeNode)\r(freeNode);
    %u(freeNode)=amg(At(freeNode,freeNode),r(freeNode));
end