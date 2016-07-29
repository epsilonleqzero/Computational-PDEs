%% Finite Element Method - Transient Maxwell Equations
%
%   Written by Ted Kwan for Math 226C Project Maxwell
%
%   Returns u(x,y,z) and p the Lagrange multiplier.
%
%   This function uses a multi-grid method combined with
%   distributed gauss-seidel to solve the equation.
%
%%% Inputs
%
% * node - Nx2 Matrix of node coordinates.
% * elem - NTx3 Matrix of elements.
% * HBo - Initial HB matrix.
% * bdFlag - Nx1 Matrix holding the boundary conditions.
% * pde.f - Function handle for the right hand side.
% * pde.g_D - Function handle for the boundary condtions.
% * J - Integer number of levels for V-cycle.
% * tol - Tolerance to solve below.
% * itrmax - max number of iterations used.
%
function [u,p,node,elem,edge,elem2edge,A,b,B,k,tol]= FiniteElemMaxwellMG(node,elem,HBo,bdFlag,pde,J,tol,itrmax)
    %% Initial Setup
    %
    HB=cell(J,1); N=zeros(J,1); NT=zeros(J,1);
    NE=zeros(J,1); HB{1}=HBo;
    [elem2edge,edge,elem2edgeSign] = dof3edge(elem);
    f=pde.f;
    %% Assemble Matrices
    %
    [A,B,M,Ap,volume,edge,locEdge,Dlambda]=AssembleStiffnessMG(node,elem);
    freeNodes=cell(J,1); BdNodes=cell(J,1); As=cell(J,1); Aps=cell(J,1);
    freeEdges=cell(J,1); BdEdges=cell(J,1); Bs=cell(J,1); Ms=cell(J,1);
    isBdEdges=cell(J,1); isBdNodes=cell(J,1);
    elems=cell(J,1); elems{1}=elem;
    edges=cell(J,1); edges{1}=edge;
    nodes=cell(J,1); nodes{1}=node;
    [bdEdge,bdNode,isBdEdges{1},isBdNodes{1}] = findboundaryedge3(edge,elem2edge,bdFlag);
    freeEdges{1}= find(~isBdEdges{1}); freeNodes{1}= find(~isBdNodes{1});
    %% First Step
    %
    N(1) = size(node,1);
    NT(1) = size(elem,1);  NE(1) = size(edge,1);
    As{1}=A(freeEdges{1},freeEdges{1}); clear A;
    Bs{1}=B(freeNodes{1},freeEdges{1}); clear B;
    Ms{1}=M(freeNodes{1},freeNodes{1}); clear M;
    Pro=cell(J,1); Res=cell(J,1); Smoothers=cell(J,4);
    Prop=cell(J,1); Resp=cell(J,1);
    %% Refine Grid and Save
    %
    for i=2:J
        [nodes{i},elems{i},bdFlag,HB{i}] = uniformrefine3(nodes{i-1},elems{i-1},bdFlag,HB{i-1});
        N(i)=size(nodes{i},1);
        [elem2edge,edges{i},elem2edgeSign] = dof3edge(elems{i});
        NT(i) = size(elems{i},1);  NE(i) = size(edges{i},1);
        [bdEdge,bdNode,isBdEdges{i},isBdNodes{i}] = findboundaryedge3(edges{i},elem2edge,bdFlag);
        freeNodes{i}=find(~isBdNodes{i});
        BdNodes{i}=bdNode;
        freeEdges{i}=find(~isBdEdges{i});
        BdEdges{i}=bdEdge;
        Protemp = transferedgered3(elems{i-1},elems{i});
        Pro{i}=Protemp(freeEdges{i},freeEdges{i-1});
        Res{i}=Pro{i}';
    end
    
    %% Construct Matrices
    %
    [HB1, NL, level] = HBstructure3(elems{J},N(1));
    [Protmp,Restmp] = transferoperator(HB1,NL);
    for i=1:J
        [A,B,M,Ap,volume,edge,LexEdge,Dlambda]=AssembleStiffnessMG(nodes{i},elems{i});
        As{i}=A(freeEdges{i},freeEdges{i});
        Ap=Ap(freeNodes{i},freeNodes{i});
        Bs{i}=B(freeNodes{i},freeEdges{i});
        Ms{i}=M(freeNodes{i},freeNodes{i});
        if(i<J)
            Proptemp=Protmp{i};
            Prop{i}=Proptemp(freeNodes{i+1},freeNodes{i});
        end
        if(i>1)
            Resptemp=Restmp{i};
            Resp{i}=Resptemp(freeNodes{i-1},freeNodes{i});
        end
        Smoothers{i,1}=tril(As{i});
        Smoothers{i,2}=triu(As{i});
        Smoothers{i,3}=tril(Ap);
        Smoothers{i,4}=triu(Ap);
    end
    b = FiniteElemMWRHS(nodes{J},elems{J},elem2edge,elem2edgeSign,f,NT(J),NE(J),LexEdge,Dlambda,volume);
    u=zeros(size(b));
    u(isBdEdges{J})=FiniteElemMWBD(nodes{J},BdEdges{J},pde.g_D);
    g = -B*u; ru=(b+((M*B)'*g)-A*u);
    rp=g(freeNodes{J});
    p=zeros(size(g));
    res=zeros(itrmax,1); tol=tol*norm(ru)
    err=2*tol; k=1;
    %% Run MG
    %
    while(err>tol && k<itrmax)
        [eu,ep] = FEMvCycleMWF(ru(freeEdges{J}),rp,J,As,Bs,Ms,Res,Pro,Prop,Resp,Smoothers,1);
        u(freeEdges{J})=u(freeEdges{J})+eu;
        p(freeNodes{J})=p(freeNodes{J})+ep;
        g = -B*u; %g=g-mean(g); 
        rp=g(freeNodes{J});
        ru=(b+((M*B)'*g)-A*u); err=norm(ru)
        res(k)=err; k=k+1;
    end
    res=res(find(res));
    node=nodes{J}; elem=elems{J}; edge=edges{J};
    %% Setup to Solve.
    %
end