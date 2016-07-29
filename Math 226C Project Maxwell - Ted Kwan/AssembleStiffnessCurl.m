%% Assemble Stiffness Matrix for the Curl-Curl Operator.
%
%   Written By Ted Kwan for Math 226C
%
%   Quick method to generate sparse matrix
%   A, the stiffness matrix, and B the divergence matrix.
%
%   This method also outputs the volume, gradient basis, edges,
%   and the lexographic ordering.
%
function [A,B,volume,LexEdge,Dlambda] = AssembleStiffnessCurl(node,elem)
    %%% Initial Setup
    %
    [elem2edge,edge,elem2edgeSign] = dof3edge(elem);
    [Dlambda,volume,~] = gradbasis3(node,elem);
    N=size(node,1); NT=size(elem,1); NE = size(edge,1);
    [Nx,Ny,~]=size(Dlambda);
    curlPhi=zeros(Nx,Ny,6); phi=zeros(NT,Ny,6);
    %%% Edge lexicographic ordering:
    %
    %   [1 2], [1 3], [1 4], [2 3], [2 4], [3 4]
    LexEdge = [1 2;...
               1 3;...
               1 4;...
               2 3;...
               2 4;...
               3 4];
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
    %% Assemble Matrices
    %
    %%% Create Divergence Matrix
    %
    index=0; indexB=0;
    for j=1:6
        for i=1:4
            
            Bii(indexB+1:indexB+NT)=double(elem(:,i));
            Bjj(indexB+1:indexB+NT)=double(elem2edge(:,j));
            sB(indexB+1:indexB+NT)=volume.*dot(phi(:,:,j),Dlambda(:,:,i),2).*elem2edgeSign(:,j);
            indexB=indexB+NT;
        end
    end
    %%% Create Stiffness Matrix
    %
    for i=1:6
       for j=1:6
           ii(index+1:index+NT)=double(elem2edge(:,i));
           jj(index+1:index+NT)=double(elem2edge(:,j));
           Aij=dot(curlPhi(:,:,i),curlPhi(:,:,j),2).*volume;
           Aij=Aij.*elem2edgeSign(:,j).*elem2edgeSign(:,i);
           sA(index+1:index+NT)=Aij;
           index=index+NT;
       end
    end
    %% Put together sparse Matrices
    %
    A=sparse(ii,jj,sA,NE,NE);
    B=sparse(Bii,Bjj,sB,N,NE);
end