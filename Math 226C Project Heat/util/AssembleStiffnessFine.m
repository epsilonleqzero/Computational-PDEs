%% Assemble Stiffness Matrix.
%
%   Quick method to generate sparse matrix
%   A, the stiffness matrix.
%
function [A,M,area] = AssembleStiffnessFine(node,elem)
    %%% Create Stiffness Matrix
    %
    N=size(node,1); NT=size(elem,1);
    ii=zeros(9*NT,1); jj=zeros(9*NT,1); sA=zeros(9*NT,1);
    ve(:,:,3)=node(elem(:,2),:)-node(elem(:,1),:);
    ve(:,:,1)=node(elem(:,3),:)-node(elem(:,2),:);
    ve(:,:,2)=node(elem(:,1),:)-node(elem(:,3),:);
    area=0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
    index=0;
    for i=1:3
       for j=1:3
           ii(index+1:index+NT)=elem(:,i);
           jj(index+1:index+NT)=elem(:,j);
           sA(index+1:index+NT)=dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
           index=index+NT;
       end
    end
    A=sparse(ii,jj,sA,N,N);
    %%% Create Mass Matrix by Mass Lumping
    %
    M = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
end