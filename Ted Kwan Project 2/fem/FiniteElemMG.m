%% Finite Element Method - Multigrid method
%
%   Written by Ted Kwan for Math 226B
%
%   Returns u(x,t) to approximate the solution
%   to the Laplace equation -\Delta u=f
%
%   This function runs the finite element method
%   using the multi-grid method to solve the linear
%   system. The function can also run the preconditioned
%   conjugate gradient method by adding pcg as a 4th
%   argument.
function [u,node,elem,k,res,tol]= FiniteElemMG(node,elem,f,J,varargin)
    %% Initial Setup
    %
    HB=cell(J,1); N=zeros(J+1,1);
    freeNodes=cell(J,1); BdNodes=cell(J,1); As=cell(J,1);
    %% Refine Grid and Save
    %
    for i=1:J
        N(i)=length(node(:,1));
        [~,~,isBdNode]=findboundary(elem);
        freeNodes{i}=find(~isBdNode);
        BdNodes{i}=find(isBdNode);
        [node,elem,bdFlag,HB{i}] = uniformrefine(node,elem);
    end
    N(J+1)=length(node(:,1));
    [bdNode,~,isBdNode]=findboundary(elem);
    freeNode=find(~isBdNode);
    [Asfull,area]=AssembleStiffnessFine(node,elem);
    As{J}=Asfull(freeNode,freeNode); clear Asfull;
    Pro=cell(J,1); Res=cell(J,1); Smoothers=cell(J+1,2);
    %% Prol, Res, Stiff And Smoother Matrices
    %
    for i=J:-1:2
        Protemp=ProHB(N(i),N(i+1),HB{i});
        if(i==J)
            Pro{i}=Protemp(freeNode,freeNodes{i});
        else
            Pro{i}=Protemp(freeNodes{i+1},freeNodes{i});
        end
        Res{i}=Pro{i}';
        As{i-1}=Res{i}*As{i}*Pro{i};
        Smoothers{i,1}=tril(As{i});
        Smoothers{i,2}=triu(As{i});
    end
    %% Calculate RHS
    %
    b=FiniteElemRHS(node,elem,f,area,N(J+1));
    %% Boundary Condition.
    %
    %%% Dirichlet Boundary Conditions.
    %
    g_D=varargin{1}; u = zeros(N(J+1),1); r = b(freeNode);
    u(bdNode) = g_D(node(bdNode,1),node(bdNode,2));
    %r(bdNode) = g_D(node(bdNode,1),node(bdNode,2));
    itrmax=varargin{2}; tol=varargin{3};
    %% Solve Au=b
    %
    if(length(varargin)>3 && strcmpi(varargin{4},'pcg'))
        %%% Preconditioned CG
        %
        [u(freeNode),k,res]=PCGvCycleAlg(r,u(freeNode),J,As,...
                            Res,Pro,Smoothers,tol,itrmax,1);
    else
        %%% Multi-Grid Method
        %
        res=zeros(itrmax,1); tol=tol*norm(r);
        err=2*tol; k=1;
        while(err>tol && k<itrmax)
        e=FEMvCycleAlg(r,J,As,Res,Pro,Smoothers,1);
        u(freeNode)=u(freeNode)+e;
        r=(r-As{J}*e); err=norm(r);
        res(k)=err; k=k+1;
        end
        res=res(find(res));
    end
end