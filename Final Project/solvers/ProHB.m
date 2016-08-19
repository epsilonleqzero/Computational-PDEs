%% Create Prolongation Matrix
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements constructs the
%   prolongation matrix using a Hierarchical basis.
%   Written using notes from Professor Chen.
function [Pro] = ProHB(nCoarseNode,nTotal,HB)
    nFineNode=nTotal-nCoarseNode; %Number of fine nodes.
    coarseNode=[1:nCoarseNode]'; coarseNodeFineIdx=coarseNode;
    ii=[coarseNodeFineIdx; double(HB(:,1)); double(HB(:,1))];
    jj=[coarseNode; double(HB(:,2)); double(HB(:,3))];
    ss=[ones(nCoarseNode,1);0.5*ones(nFineNode, ...
        1);0.5*ones(nFineNode,1)];
    Pro=sparse(ii,jj,ss,nTotal,nCoarseNode);
end