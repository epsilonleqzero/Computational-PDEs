%% Set the boundary flag
%
%   Written by Ted Kwan for Math 226B  Project 1.
%
%   This method sets the boundary flag for Neumann
%   Boundary conditions. It checks each boundary element
%   then sets the non-boundary node to 2 so that the
%   opposite edge will be counted as a boundary edge.
%
function [bdFlag] = SetBdFlagCir(isBdElem,elem,isBdNode,NT,d)
    bdFlag=int8(zeros(NT,d+1));
    for i=1:NT
       if(isBdElem(i))
        %%% Find Boundary node
        %
        %   Once the node is found, then it is checked to see if there
        %   are two boundary nodes. If there are, the flag is set.
        if(((isBdNode(elem(i,1))==1) && (isBdNode(elem(i,2))==1)) ...
            || ((isBdNode(elem(i,1))==1) && (isBdNode(elem(i,3))==1))...
            || ((isBdNode(elem(i,2))==1) && (isBdNode(elem(i,3))==1)))
            %%% Change bdFlag for a boundary node.
            for j=1:(d+1)
                if(isBdNode(elem(i,j))==0)
                    bdFlag(i,j)=2;
                end
            end
        end
       end
    end
end