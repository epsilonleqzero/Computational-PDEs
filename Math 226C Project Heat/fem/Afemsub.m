%% Adaptive Finite Element Subroutine
%
%   Written by Ted Kwan for Math 226C Project 1.
%
%   This function performs the adaptive finite element
%   on for each time step. It performs a different calculation
%   based on the residual used.
%
%   This function does both refinement and coarsening on the
%   mesh at each time step.

function [u,M,A,At,node,elem,err,h1err,Ns] = Afemsub(node,elem,pde,uold,dt,icflag,varargin)
    %% Initial Setup
    %
    type=varargin{1};
    if(type==1)
        maxN = 5e3; theta = 0.6; theta2 = 0.08*theta; maxIt = 7;
    elseif(type==2)
        maxN = 5e3; theta = 0.4; theta2 = 0.2*theta ; maxIt = 8;
    elseif(type==3)
        maxN = 5e3; theta = 0.6; theta2 = 0.1*theta ; maxIt = 8;
    else
        maxN = 5e3; theta = 0.6; theta2 = 0.3*theta ; maxIt = 8;
    end
    N = zeros(maxIt,1); Ns = zeros(maxIt,1);
    g_D=@pde.g_D; f=@pde.f; f1=@pde.f1;
    err=zeros(maxIt,1); h1err=zeros(maxIt,1);
    err(1)=1; N(1)=size(node,1);
    %figure(1); 
    %%  Adaptive Finite Element Method
    %
    % *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE/COARSEN*
    for k=1:maxIt
        %% Step 1: SOLVE
        %
        [u,M,A,At]=afemheatsolve(node,elem,g_D,f1,uold,dt);
        %% Plot mesh and solution
        %showresult(node,elem,u,[-9,15]);
        %axis([0,1,0,1,-0.01,0.1]); pause(0.001);
        %showmesh(node,elem); pause(0.0001);
        %% Step 2: Estimate Residual
        %
        pde1=pde;
        pde1.f=@(x) (f(x(:,1),x(:,2)));
        if(type==1)
            eta = estimateresidualheat(node,elem,u,pde1,uold,dt);
        elseif(type==2)
            eta = estimateresidualheata(node,elem,u,pde1,uold,dt);
        elseif(type==3)
            eta = estimateresidualheatb(node,elem,u,pde1,uold,dt);
        else
            eta = estimaterecovery(node,elem,u);
        end
        N(k) = size(node,1); Ns(k) = length(node(:,1));
        %% Record Error
        %
        %   This is not used in the actual algorithm.
        utestL2=g_D(node(:,1),node(:,2));
        err(k)=sqrt((((utestL2-u))')*M*(utestL2-u));
        h1err(k)=sqrt((((utestL2-u))')*A*(utestL2-u));
        if (k==maxIt)
            break;
        end
        %% Mark and refine or coarsen
        %
        if(N(k)<maxN)
            %% Refine
            %
            markedElem = mark(elem,eta,theta);
            [node,elem,~,HB] = bisect(node,elem,markedElem);
            if(icflag>0)
                u0=varargin{2};
                uold=u0(node(:,1),node(:,2));
            else
                uold=nodeinterpolate(uold,HB);
            end
        else
            %% Coarsen
            %
            markedElem = mark(elem,eta,theta2,'COARSEN');
            [node,elem,~,indexMap,~] = coarsen(node,elem,markedElem);
            uold=nodeinterpolate(uold,indexMap);
        end
    end
end