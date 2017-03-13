%% Generate Mesh
%
%   Written by Ted Kwan for Math 226B Project 1
%
%   This method generates the meshes, plots them,
%   then compares them in the profile viewer.
%
clear all;
[node,elem] = squaremesh([0,1,0,1],0.25);
% [node,elem] = circlemesh(0,0,1,0.2);
% showmesh(node,elem);
profile on
for i=1:3
    [node,elem] = uniformrefine(node,elem);
    tic; assemblingstandard(node,elem); toc;
    tic; assemblingsparse(node,elem); toc;
    tic; assembling(node,elem); toc;
end
profile viewer