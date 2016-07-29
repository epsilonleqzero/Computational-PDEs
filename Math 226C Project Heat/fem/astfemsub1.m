%% Adaptive Finite Element Subroutine
%
%   Written by Ted Kwan for Math 226C Project 1.
%
%   This function performs the adaptive finite element
%   on for each time step. It performs a different calculation
%   based on the residual used.
%
%   This function does both refinement and coarsening on the
%   mesh at each time step. In addition, there is time adaptivity
%   when needed. This is buggy, and needs to have some issues
%   worked out.

function [u,M,A,At,node,elem,t,err,h1err,Ns] = astfemsub1(node,elem,pde,uold,dt,told,tf,icflag,varargin)
    maxN = 6e3; theta = 0.7; maxIt = 8;
    N = zeros(maxIt,1); Ns = zeros(maxIt,1);
    g_D=@pde.g_D; f=@pde.f;
    err=zeros(maxIt,1); h1err=zeros(maxIt,1);
    err(1)=1; N(1)=size(node,1); stop=1;
    %figure(1); 
    %%  Adaptive Finite Element Method
    % *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE/COARSEN*
    
    for k=1:maxIt
        %% Step 1: SOLVE
        %
        pde1.f=@(x) (f(x(:,1),x(:,2),told+(dt/2)));
        pde1.g_D=@(x,y) pde.g_D(x,y,told+dt);
        ft=@(x,y) (f(x,y,told)+f(x,y,told+dt))./2;
        [u,M,A,At]=afemheatsolve(node,elem,pde1.g_D,ft,uold,dt);        
        %% Plot mesh and solution
        %
        %figure(1);  showresult(node,elem,u,[-9,15]);
        %axis([0,1,0,1,-0.01,0.1]);
        showmesh(node,elem); pause(0.005);
        %% Step 2: ESTIMATE
        %
        eta = astestimateresidualheat(node,elem,u,pde1,uold,dt);
        N(k) = size(node,1); Ns(k) = length(node(:,1));
        utestL2=g_D(node(:,1),node(:,2),told+dt);
        err(k)=sqrt((((utestL2-u))')*M*(utestL2-u));
        h1err(k)=sqrt((((utestL2-u))')*A*(utestL2-u));
        if (k==maxIt)
            break;
        end
        tol=1e-3; checkt=norm(eta);
        if(checkt>10*tol && stop<5)
            stop=stop+1
            dt=dt/2
        else
            if(N(k)<maxN || (checkt>tol && N(k)<1.5*maxN))
                %% Mark
                %
                markedElem = mark(elem,eta,theta);
                %% Refine
                %
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
                markedElem = mark(elem,eta,0.2*theta,'COARSEN');
                [node,elem,~,indexMap,~] = coarsen(node,elem,markedElem);
                uold=nodeinterpolate(uold,indexMap); 
            end
        end
    end
    t=told+dt;
end