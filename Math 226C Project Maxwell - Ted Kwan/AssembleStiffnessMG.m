%% Assemble Stiffness Matrix for the DGS and MG methods.
%
%   Written By Ted Kwan for Math 226C
%
%   Quick method to generate sparse matrix
%   A, the stiffness matrix for the discrete Laplacian, and 
%   B the divergence matrix. It also creates the matrix for the
%   discrete poisson problem with homogeneous Neumann boundary 
%   conditions.
%
%   This method also outputs the volume, gradient basis, edges,
%   and the lexographic ordering.
%
function [A,B,M,Ap,volume,edge,LexEdge,Dlambda] = AssembleStiffnessMG(node,elem)
    %%% Initial Setup
    %
    [elem2edge,edge,elem2edgeSign] = dof3edge(elem);
    [Dlambda,volume,~] = gradbasis3(node,elem);
    N=size(node,1); NT=size(elem,1); NE = size(edge,1);
    %%% Edge lexicographic ordering:
    %
    %   [1 2], [1 3], [1 4], [2 3], [2 4], [3 4]
    LexEdge = [1 2;...
               1 3;...
               1 4;...
               2 3;...
               2 4;...
               3 4];
    [Nx,Ny,~]=size(Dlambda);
    curlPhi=zeros(Nx,Ny,6); phi=zeros(NT,Ny,6);
    %% Calculate Curl(Phi)
    %
    %   curl(\phi_k)=2\nabla lambda_i \times \nabla lambda_j
    %
    for i=1:6
        li = LexEdge(i,1); lj = LexEdge(i,2);
        curlPhi(:,:,i)=2*cross(Dlambda(:,:,li),Dlambda(:,:,lj),2);
        phi(:,:,i) = ((Dlambda(:,:,lj)-Dlambda(:,:,li)))/4;
    end
    Bii=zeros(24*NT,1); Bjj=zeros(24*NT,1); sB=zeros(24*NT,1);
    sA=zeros(36*NT,1); ii=zeros(36*NT,1); jj=zeros(36*NT,1);
    %% Assemble Stiffness Matrix
    %
    index=0; indexB=0;
    for i=1:6
       for j=1:6
           ii(index+1:index+NT)=double(elem2edge(:,i));
           jj(index+1:index+NT)=double(elem2edge(:,j));
           Aij=dot(phi(:,:,i),phi(:,:,j),2).*volume;
           Aij=Aij.*elem2edgeSign(:,j).*elem2edgeSign(:,i);
           sA(index+1:index+NT)=Aij;
           index=index+NT;
       end
    end
    %% Create Mass Matrix by Mass Lumping
    %
    Md = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],[volume;volume;volume;volume]/4,[N,1]);
    %Md=Md.^-1;
    %% Create Divergence Matrix
    %
    for j=1:6
        for i=1:4
            
            Bii(indexB+1:indexB+NT)=double(elem(:,i));
            Bjj(indexB+1:indexB+NT)=double(elem2edge(:,j));
            sB(indexB+1:indexB+NT)=volume.*dot(phi(:,:,j),Dlambda(:,:,i),2).*elem2edgeSign(:,j);
            indexB=indexB+NT;
        end
    end
    %% Create Stiffness Matrix for Poisson Problem in DGS
    %
    sAp=zeros(10*NT,1); Pii=zeros(10*NT,1); Pjj=zeros(10*NT,1);
    indexp=0;
    for i = 1:4
        for j = i:4
           Pii(indexp+1:indexp+NT)=(elem(:,i));
           Pjj(indexp+1:indexp+NT)=(elem(:,j));
           Aij=dot(Dlambda(:,:,i),Dlambda(:,:,j),2).*volume;
           sAp(indexp+1:indexp+NT)=Aij;
           indexp=indexp+NT;      
        end
    end
    Ap=sparse(Pii,Pjj,sAp,N,N);
    %%% Fixup for Neumann BC
    %
    fixedNode = 1;
    bdidx = zeros(N,1); 
    bdidx(fixedNode) = 1;
    Tbd = spdiags(bdidx,0,N,N);
    T = spdiags(1-bdidx,0,N,N);
    Ap = T*Ap*T + Tbd;
    %% Sparse Matrix Assembling
    %
    A=sparse(ii,jj,sA,NE,NE);
    M=spdiags(Md,0,N,N);
    B=sparse(Bii,Bjj,sB,N,NE);
end