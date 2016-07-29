%% Finite Element Method BD
%
%   Written by Ted Kwan for Math 226C Project Maxwell
%
%   This function calculates the boundary on edge
%   elements for the finite element method.
%
function [u] = FiniteElemMWBD(node,edge,g_D)
    %% Calculate Boundary Values using Simpson's method.
    %
    ge1 = g_D(node(edge(:,1),1),node(edge(:,1),2),node(edge(:,1),3));
    ge2 = g_D(node(edge(:,2),1),node(edge(:,2),2),node(edge(:,2),3));
    mid=(node(edge(:,1),:)+node(edge(:,2),:))/2;
    gmid = g_D(mid(:,1),mid(:,2),mid(:,3));
    edgev=node(edge(:,2),:)-node(edge(:,1),:);
    %%% Simpson's Method
    %
    u=dot(edgev,(ge1+ge2+4*gmid)/6,2);
end