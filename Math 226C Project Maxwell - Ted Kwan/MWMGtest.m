%% Test For Maxwell Finite Element Multigrid Method
%
%   Written by Ted Kwan for Math 226C Project Maxwell.
%
%   This script runs the maxwell equation on a cube grid
%   which is missing a section. It can be run against 3
%   different test functions.
%
%   This tests the standard finite element discretization, 
%   using a direct-solver.
%
clear all; %close all;
%%% For DGS
%
%[nodes,elems,HB] = cubemesh([-1,1,-1,1,-1,1],0.25);

%%% For MG
%
[nodes,elems,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);

[nodes,elems] = delmesh(nodes,elems,'x<0 & y<0 & z>0');

%% Test 1
%
pde.f=@(x,y,z) [sin(z),cos(x),sin(y)];
pde.g_D=@(x,y,z) [sin(z),cos(x),sin(y)];
pde.g_Dmg=@(x,y,z) [zeros(size(x(:,1))),zeros(size(x(:,2))),zeros(size(x(:,3)))];
pde.g_D1=@(x) [sin(x(:,3)),cos(x(:,1)),sin(x(:,2))];
pde.curlu=@(x) [cos(x(:,2)),cos(x(:,3)),-sin(x(:,1))];

% %% Test 2
% %
% pde.f=@(x,y,z) -2*[ones(size(z)),ones(size(x)),ones(size(y))];
% pde.g_D=@(x,y,z) [z.^2,x.^2,y.^2];
% pde.g_D1=@(x) [(x(:,3).^2),x(:,1).^2,x(:,2).^2];
% pde.curlu=@(x) 2*[(x(:,2)),((x(:,3))),((x(:,1)))];

% %% Test 3
% %
% pde.g_D=@(x,y,z) [0*x,cos(x), cos(x)];
% pde.f=@(x,y,z) [0*x,cos(x),cos(x)];
% pde.g_D1=@(x) [0*(x(:,1)),cos(x(:,1)),cos(x(:,1))];
% pde.curlu=@(x) [0*x(:,1), sin(x(:,1)), -sin(x(:,1))];

%% Start tests
%
bdFlag = setboundary3(nodes,elems,'Dirichlet');
maxIt=3;
Js=[2,3,4]';
hs=[0.25,0.125,0.0625]';
L2Err = zeros(maxIt,1); L2Err2 = zeros(maxIt,1);
N = zeros(maxIt,1); uIuhErr = zeros(maxIt,1);
pErr = zeros(maxIt,1); pErr1 = zeros(maxIt,1);
itrs = zeros(maxIt,1); tols = zeros(maxIt,1);
for i=1:maxIt
    [ut,p,node,elem,edge,elem2edge,A,b,B,itrs(i),tols(i)]=FiniteElemMaxwellMG(nodes,elems,HB,bdFlag,pde,Js(i),hs(i)^2,50);
    %[ut,p,node,elem,edge,elem2edge,A,b,B,itrs(i),tols(i)]=FiniteElemMaxwellDGS(nodes,elems,bdFlag,pde,hs(i)^2,1000);
     uI = edgeinterpolate(pde.g_D1,node,edge);
     L2Err(i) = getL2error3NE(node,elem,pde.g_D1,(ut));
     pErr(i)=(norm(B*ut));
     pErr1(i)=(norm(B'*p));
     uIuhErr(i) = sqrt((ut-uI)'*(A)*(ut-uI));
     N(i) = length(ut);
     %%% Uncomment for DGS
     %
     %[nodes,elems,bdFlag] = uniformrefine3(nodes,elems,bdFlag);
end
mgt=[hs,Js,tols,itrs,uIuhErr,L2Err];
dlmwrite('mgtests.txt', mgt , '&');
figure;
showrate2(N,L2Err,2,'b-+','||u-u_h||',...
          N,uIuhErr,2,'r-+','||uI-u_h||_A');
figure;
showrate2(N,pErr,2,'b-+','||\nabla \cdot u_h||',...
          N,pErr1,2,'r-+','||\nabla \cdot p||');