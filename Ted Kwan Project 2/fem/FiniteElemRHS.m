%% Finite Element Method RHS
%
%   Written by Ted Kwan for Math 226B Project 2
%   This function calculates the right-hand side of the FEM.
function [b] = FiniteElemRHS(node,elem,f,area,N)
    %% Calculate RHS
    mid1 = (node(elem(:,2),:)+node(elem(:,3),:))/2;
    mid2 = (node(elem(:,3),:)+node(elem(:,1),:))/2;
    mid3 = (node(elem(:,1),:)+node(elem(:,2),:))/2;
    bt1 = area.*(f(mid2(:,1),mid2(:,2))+f(mid3(:,1),mid3(:,2)))/6;
    bt2 = area.*(f(mid3(:,1),mid3(:,2))+f(mid1(:,1),mid1(:,2)))/6;
    bt3 = area.*(f(mid1(:,1),mid1(:,2))+f(mid2(:,1),mid2(:,2)))/6;
    b = accumarray(elem(:),[bt1;bt2;bt3],[N,1]);
end

