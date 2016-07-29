%% Lid Driven Cavity Boundary Conditions
%
%   Written by Ted Kwan for Math 226C Project 2.
%
%   This function calculates the boundary conditions
%   for the lid driven cavity problem.
%
function [u] = LidCavityu(x,y)
    if(y==1)
        u=1;
    else
        u=0;
    end
end

