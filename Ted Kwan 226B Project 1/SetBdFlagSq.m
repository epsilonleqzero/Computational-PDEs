%% Set the boundary flag
%
%   Written by Ted Kwan for Math 226B  Project 1.
%
%   This method sets the boundary flag for Neumann
%   Boundary conditions. It checks each boundary element
%   then sets the node opposite each edge to 2.
%
function [bdFlag] = SetBdFlagSq(isBdElem,elem,isBdNode,bdEdge,NT,d)
    bdFlag=int8(zeros(NT,d+1));
    Nbde=length(bdEdge(:,1));
    for i=1:NT
       if(isBdElem(i))
            %%% Boundary Element has been found.
            %
            %   checks to ensure that the edges are
            %   actually on the boundary, then sets node.
            el=elem(i,:);
            if(ismember([el(1),el(2)],bdEdge,'rows') ||...
               ismember([el(2),el(1)],bdEdge,'rows'))
               bdFlag(i,3)=2;
            end
            if(ismember([el(1),el(3)],bdEdge,'rows') ||...
                   ismember([el(3),el(1)],bdEdge,'rows'))
                bdFlag(i,2)=2;
            end
            if(ismember([el(2),el(3)],bdEdge,'rows') ||...
                   ismember([el(3),el(2)],bdEdge,'rows'))
                bdFlag(i,1)=2;
            end
       end
    end
end
