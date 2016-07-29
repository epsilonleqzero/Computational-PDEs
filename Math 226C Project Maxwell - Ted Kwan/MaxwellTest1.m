%% Test For Maxwell Finite Element Method
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
%% Mesh Setup
%
[node,elem] = cubemesh([-1,1,-1,1,-1,1],0.5);
[node,elem] = delmesh(node,elem,'x<0 & y<0 & z>0');
bdFlag = setboundary3(node,elem,'Dirichlet');
maxIt=3;

% %% Test 1
% %
% pde.f=@(x,y,z) [sin(z),cos(x),sin(y)];
% pde.g_D=@(x,y,z) [sin(z),cos(x),sin(y)];
% pde.g_D1=@(x) [sin(x(:,3)),cos(x(:,1)),sin(x(:,2))];
% pde.curlu=@(x) [cos(x(:,2)),cos(x(:,3)),-sin(x(:,1))];

%% Test 1
%
pde.f=@(x,y,z) [zeros(size(z)),-5*ones(size(x)),zeros(size(y))];
pde.g_D=@(x,y,z) [x.*y,2*x.*z,6*y.*z];
pde.g_D1=@(x) [x(:,1).*x(:,2),2*x(:,1).*x(:,3)...
               ,6*x(:,2).*x(:,3)];
pde.curlu=@(x) [6*x(:,3)-2*x(:,1),zeros(size(x(:,1))),...
                2*x(:,3)-x(:,1)];
pde.divu=@(x) 8*x(:,2);
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

%% Error Recording Setup.
%
L2Err = zeros(maxIt,1); L2Err2 = zeros(maxIt,1);
N = zeros(maxIt,1); uIuhErr = zeros(maxIt,1);
pErr = zeros(maxIt,1);pErr1 = zeros(maxIt,1);
for i=1:maxIt
    %% Main Loop
    %
    [ut,p,node,elem,edge,elem2edge,A,b,B]=FiniteElemMaxwell(node,elem,bdFlag,pde);
    uI = edgeinterpolate(pde.g_D1,node,edge);
    L2Err(i) = getL2error3NE(node,elem,pde.g_D1,(ut));
    divu=pde.divu(node);
    pErr(i)=(norm(divu-B*ut));
    %pErr(i)=sqrt((ut-uI)'*(B)*(ut-uI));
    pErr1(i)=(norm(B'*p));
    %uIuhErr(i) = sqrt((ut-uI)'*(A)*(ut-uI));
    uIuhErr(i) = sqrt((ut-uI)'*(A)*(ut-uI));
    %uIuhErr(i) = getHcurlerror3NE(node,elem,pde.curlu,real(ut));
    N(i) = length(ut);
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
end
%% Plot Results
%
%%% L2 Error and H1 error for u
%
figure;
showrate2(N,L2Err,2,'b-+','||u-u_h||',...
          N,uIuhErr,2,'r-+','||uI-u_h||_A');
%%% Divergence error
%
figure;
showrate2(N,pErr,2,'b-+','||\nabla \cdot u_h||',...
          N,pErr1,2,'r-+','||\nabla \cdot p||');