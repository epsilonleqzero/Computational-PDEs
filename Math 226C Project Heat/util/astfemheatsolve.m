%% AFEM Heat Solving method
%
%
function [u,Md,A,At] = astfemheatsolve(node,elem,g_D,f,uold,dt,told)
    [A,M,area]=assemblematrix(node,elem);
    [bdNode,~,isBdNode]=findboundary(elem);
    freeNode= find(~isBdNode);
    N=length(node(:,1)); Md=spdiags(M,0,N,N);
    At=(Md+(dt/2).*A);
    Atr=((Md)-(dt/2).*A);
    u=zeros(N,1);
    u(bdNode) = g_D(node(bdNode,1),node(bdNode,2));
    br=FiniteElemRHS(node,elem,f,area,N);
    %%% Create RHS
    %
    br=(dt).*br; b1=Atr*uold;
    b=b1+br; r=b-At*u;
    %%% End RHS
    u=amg()
    u(freeNode)=At(freeNode,freeNode)\r(freeNode);
end