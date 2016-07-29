%% Finite Element Method - Transient Maxwell Equations
%
%   Written by Ted Kwan for Math 226C Project Maxwell
%
%   Returns u(x,y,z) and p the Lagrange multiplier.
%
%   This function uses the distributed gauss-seidel method
%   to approximate the solution.
%
%%% Inputs
%
% * node - Nx2 Matrix of node coordinates.
% * elem - NTx3 Matrix of elements.
% * bdFlag - Nx1 Matrix holding the boundary conditions.
% * pde.f - Function handle for the right hand side.
% * pde.g_D - Function handle for the boundary condtions.
% * tol - Tolerance to solve below.
% * itrmax - max number of iterations used.
%
function [u,p,node,elem,edge,elem2edge,A,b,B,k,tol]= FiniteElemMaxwellDGS(node,elem,bdFlag,pde,tol,itrmax)
    %% Initial Setup
    %
    [elem2edge,edge,elem2edgeSign] = dof3edge(elem);
    N = size(node,1);   NT = size(elem,1);  NE = size(edge,1);
    f=pde.f;
    %% BC
    %
    [bdEdge,bdNode,isBdEdge,isBdNode] = findboundaryedge3(edge,elem2edge,bdFlag);
    freeEdge= find(~isBdEdge); freeNode= find(~isBdNode);
    
    %% Assemble Matrices
    %
    [A,B,M,Ap,volume,edge,locEdge,Dlambda]=AssembleStiffnessMG(node,elem);
    b = FiniteElemMWRHS(node,elem,elem2edge,elem2edgeSign,f,NT,NE,locEdge,Dlambda,volume);
    %% Setup to Solve.
    %
    u=zeros(size(b));
    Af=A(freeEdge,freeEdge); Bf  = B(freeNode,freeEdge);
    Apf=Ap(freeNode,freeNode);
    R=tril(Af);
    Rp=tril(Apf);
    u(isBdEdge)=FiniteElemMWBD(node,bdEdge,pde.g_D);
    g = -B*u; rp=g; p=zeros(size(g));
    ru=(b+(M*B)'*g-A*u);
    res=zeros(itrmax,1); tol=tol*norm(ru)
    err=2*tol; k=1;
    %eh = triu(A)\(D.*(tril(A)\r));
    while(err>tol && k<itrmax)
        eu=zeros(size(ru(freeEdge))); ep=zeros(size(rp(freeNode)));
        [eu,ep] = DGSMaxwell(ru(freeEdge),rp(freeNode),ep,Bf,R,Rp);
        u(freeEdge)=u(freeEdge)+eu;
        g = -B*u; rp=g;
        p(freeNode)=p(freeNode)+ep;
        ru=(b+((M*B)'*g)-B'*p-A*u);
        err=norm(ru)
        res(k)=err; k=k+1;
    end
    res=res(find(res));
end