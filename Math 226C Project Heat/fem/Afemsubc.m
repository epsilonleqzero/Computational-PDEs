%% Adaptive Finite Element Subroutine
%
%   Written by Ted Kwan for Math 226C Project 1.
%
%   This function performs the adaptive finite element
%   on for each time step.
%
%   This function does only refinement on the mesh at
%   each time step until the max number of nodes has
%   been reached.

function [u,M,A,At,node,elem,err,h1err,Ns] = Afemsubc(node,elem,pde,uold,dt,icflag,varargin)
    %% Initial Setup
    %
    maxN = 7e3; theta = 0.4; maxIt = 6;
    N = zeros(maxIt,1); Ns = zeros(maxIt,1);
    g_D=@pde.g_D; f=@pde.f; f1=@pde.f1;
    err=zeros(maxIt,1); h1err=zeros(maxIt,1);
    err(1)=1; N(1)=size(node,1);
    %figure(1);
    if(N(1)<maxN)
    %%  Adaptive Finite Element Method
    % *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE/COARSEN*
    for k=1:maxIt
        %% Step 1: SOLVE
        %
        [u,M,A,At]=afemheatsolve(node,elem,g_D,f1,uold,dt);
        %% Plot mesh and solution
        %
        %figure(1);
        %showresult(node,elem,u,[-9,15]);
        %axis([0,1,0,1,-0.01,0.1]); pause(0.001);
        %showmesh(node,elem); pause(0.0001);
        %% Step 2: ESTIMATE
        %
        pde1=pde;
        pde1.f=@(x) (f(x(:,1),x(:,2)));
        eta = estimateresidualheat(node,elem,u,pde1,uold,dt);
        N(k) = size(node,1); Ns(k) = length(node(:,1));
        utestL2=g_D(node(:,1),node(:,2));
        err(k)=sqrt((((utestL2-u))')*M*(utestL2-u));
        h1err(k)=sqrt((((utestL2-u))')*A*(utestL2-u));
        if (k==maxIt)
            break;
        end
        if(N(k)<maxN)
            %% Step 3: MARK
            %
            markedElem = mark(elem,eta,theta);
            %% Step 4: REFINE
            %
            [node,elem,~,HB,~] = bisect(node,elem,markedElem);
            if(icflag>0)
                u0=varargin{2};
                uold=u0(node(:,1),node(:,2));
            else
                uold=nodeinterpolate(uold,HB);            
            end
        end
    end
    else
        [u,M,A,At]=afemheatsolve(node,elem,g_D,f1,uold,dt);
        Ns=ones(maxIt,1)*N(1); utestL2=g_D(node(:,1),node(:,2));
        err=ones(maxIt,1)*sqrt((((utestL2-u))')*M*(utestL2-u));
        h1err=ones(maxIt,1)*sqrt((((utestL2-u))')*A*(utestL2-u));
    end
end