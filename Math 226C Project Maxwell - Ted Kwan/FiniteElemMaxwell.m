%% Finite Element Method - Transient Maxwell Equations
%
%   Written by Ted Kwan for Math 226C Project Maxwell
%
%   Returns u(x,y,z) and p the Lagrange multiplier.
%
%%% Inputs
%
% * node - Nx2 Matrix of node coordinates.
% * elem - NTx3 Matrix of elements.
% * bdFlag - Nx1 Matrix holding the boundary conditions.
% * pde.f - Function handle for the right hand side.
% * pde.g_D - Function handle for the boundary condtions.
%
function [u,p,node,elem,edge,elem2edge,A,b,B]= FiniteElemMaxwell(node,elem,bdFlag,pde)
    %% Initial Setup
    %
    [elem2edge,edge,elem2edgeSign] = dof3edge(elem);
    N = size(node,1);   NT = size(elem,1);  NE = size(edge,1);
    f=pde.f;
    %% Find Boundary
    %
    [bdEdge,bdNode,isBdEdge,isBdNode] = findboundaryedge3(edge,elem2edge,bdFlag);
    freeEdge= find(~isBdEdge); freeNode= find(~isBdNode);
    
    %% Assemble Matrices
    %
    [A,B,volume,locEdge,Dlambda]=AssembleStiffnessCurl(node,elem);
    b = FiniteElemMWRHS(node,elem,elem2edge,elem2edgeSign,f,NT,NE,locEdge,Dlambda,volume);
    %% Setup to Solve.
    %
    %%% Boundary Conditions
    %
    u=zeros(size(b)); Bt=B';
    u(isBdEdge)=FiniteElemMWBD(node,bdEdge,pde.g_D);
    b=b-A*u;
    %%% Matrix Prep.
    %
    Af=A(freeEdge,freeEdge); Bf  = Bt(freeEdge,freeNode);
    Nint = size(freeNode,1); NEint = size(freeEdge,1);
    %%% Pressure Residual
    %
    g = -B*u;  p=zeros(size(g));
    g=g-mean(g);
    %% Solve Big Matrix Equation
    %
    Asys = [Af,Bf;...
            Bf',sparse(Nint,Nint)];
    tmpu = Asys\[b(freeEdge); g(freeNode)];
    u(freeEdge) = tmpu(1:NEint);
    p(freeNode) = tmpu(NEint+1:NEint+Nint);
end