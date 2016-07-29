%% Finite Element Method RHS for Maxwell Equations
%
%   Written by Ted Kwan for Math 226C Project Maxwell
%
%   This function calculates the right-hand side of 
%   finite element method for Maxwell equations on
%   edge elements.
%
function [b] = FiniteElemMWRHS(node,elem,elem2edge,elem2edgeSign,f,NT,NE,...
                                locEdge,Dlambda,volume)
    %% Calculate RHS
    %
    bt=zeros(NT,6);
    mid = (node(elem(:,1),:) ...
        + node(elem(:,2),:) ...
        + node(elem(:,3),:) ...
        + node(elem(:,4),:))/4;
    fmp=f(mid(:,1),mid(:,2),mid(:,3));
    for m=1:6
        i = locEdge(m,1); j = locEdge(m,2);
        phi = (Dlambda(:,:,j)-Dlambda(:,:,i))/4;
        bt(:,m)=dot(phi,fmp,2);
    end
    %%% Correct Sign
    %
    bt=elem2edgeSign.*bt;
    %%% One-Point Integration
    %
    bt=bt.*repmat(volume,1,6);
    b = accumarray(elem2edge(:),bt(:),[NE,1]);
end