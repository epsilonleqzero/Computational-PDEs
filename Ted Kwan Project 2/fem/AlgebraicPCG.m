%% Finite Element Method - Multigrid method
%
%   Written by Ted Kwan for Math 226B
%
%   This method solves Ax=b using an algebraic
%   preconditioned conjugate gradient method with
%   a J-level v-cycle as the preconditioning matrix.
%   The finite element method can also be solved by
%   using additional inputs.
%
%%% Inputs
%
% * node - Nx2 Matrix of node coordinates.
% * elem - NTx3 Matrix of elements.
% * J - Number of steps for v-cycle.
% * itrmax - max number of iterations.
% * tol - tolerance for solving system.
%
function [u,node,elem,stop]= AlgebraicPCG(node,elem,J,varargin)
    %% Initial Setup
    %
    [~,~,isBdNode]=findboundary(elem);
    N=length(node(:,1)); bdNode=find(isBdNode);
    freeNode= find(~isBdNode);
    Pro=cell(J,1); Res=cell(J,1); A=cell(J,1);
    itrmax=varargin{1}; Smoothers=cell(J,2);
    %% Finest Stiffness Matrix
    %
    [Afull,~,area]=assemblematrix(node,elem);
    A{J}=Afull(freeNode,freeNode);
    %% Prolongation, Restriction, Smoother and Stiffness Matrices
    %
    for i=J:-1:2
        [isC,As] = coarsenAMGc(A{i});
        [Pro{i},Res{i}] = interpolationAMGs(As,isC);
        A{i-1}=Res{i}*A{i}*Pro{i};
        Smoothers{i,1}=tril(A{i});
        Smoothers{i,2}=triu(A{i});
    end
    u = zeros(N,1);
    %% Calculate RHS and Solve
    %
    if(length(varargin)>3)
        %% Finite Element Method
        %
        f=varargin{4};
        r=FiniteElemRHS(node,elem,f,area,N);
        tol=(varargin{2});
        %%% Dirichlet Boundary Condition.
        %
        g_D=varargin{3};
        u(bdNode) = g_D(node(bdNode,1),node(bdNode,2));
        b=r(freeNode);
        %%% Solve
        %
        [u(freeNode),stop,res]=PCGvCycleAlg(b,u(freeNode),J,A, ...
                            Res,Pro,Smoothers,tol,itrmax,5);
    else
        %% Algebraic Multi-Grid
        %
        r=varargin{3}; b = r(freeNode);
        tol=varargin{2};%*norm(r);
        [u,stop,res]=PCGvCycleAlg(b,u(freeNode),J,A, ...
                            Res,Pro,Smoothers,tol,itrmax,5);
    end
end