%% Adaptive Finite Element Method - Heat Equation
%
%   Written by Ted Kwan for Math 226C
%
%   Returns u(x,t) to approximate the solution
%   to the Heat equation u_{t}-\Delta u=f Using an
%   adaptive Crank-Nicolson scheme.
%
%%% Inputs
%
% * node - Nx2 Matrix of node coordinates.
% * elem - NTx3 Matrix of elements.
% * f - Function handle for the right hand side.
% * u0 - Function handle for the initial condition.
% * t0 - Initial time.
% * tf - Final time.
% * dt - Time step size.
% * varargin options:
% 1 function handle for g(x,t).
% 2 Type of adaptive refinement to use.
% 3 Algorithm to use.
%
function [u,node,elem,t,M,A,At,errs]= afemheat(node,elem,f,u0,t0,tf,dt,varargin)
    %% Initial Setup
    %
    t=[t0:dt:tf]'; Nt=length(t);
    u=u0(node(:,1),node(:,2));
    g_D=varargin{1}; type=varargin{2};
    alg=varargin{3};
    errs=cell(Nt,4); j=1;
%     %% Capture Movie
%     % 
%     mv = VideoWriter('HeatAdaptiveAlg2.avi');
%     mv.FrameRate=3; open(mv);
%     fig=figure;
%     showresult(node,elem,u,[-9,20]); axis([0,1,0,1,-0.01,0.1]);
%     legend(['t=' num2str(t(1))],'LOCATION','best');
%     pause(0.01);
    for i=2:Nt
        uold=u;
        if(type==2)
            pde.f=@(x,y) (f(x,y,t(i))+f(x,y,t(i-1)))./2;
        else
           pde.f=@(x,y) f(x,y,t(i)-(dt/2)); 
        end
        pde.f1=@(x,y) (f(x,y,t(i))+f(x,y,t(i-1)))./2;
        pde.g_D=@ (x,y) g_D(x,y,t(i));
        if(i>=3)
            switch alg
                case 2
                    [u,M,A,At,node,elem,err,h1err,ns]= Afemsubc(node,elem,pde,uold,dt,0,type);
                otherwise
                    [u,M,A,At,node,elem,err,h1err,ns]= Afemsub(node,elem,pde,uold,dt,0,type);
            end
        else
            switch alg
                case 2
                    [u,M,A,At,node,elem,err,h1err,ns]= Afemsubc(node,elem,pde,uold,dt,1,type,u0);
                otherwise
                    [u,M,A,At,node,elem,err,h1err,ns]= Afemsub(node,elem,pde,uold,dt,1,type,u0);
            end
        end
%         %%% Uncomment to save movie
%         if(mod(i,4)==0)
%             showresult(node,elem,u,[-9,20]); axis([0,1,0,1,-0.01,0.1]);
%             legend(['t=' num2str(t(i))],'LOCATION','best');
%             pause(0.01);
%         end
%         frame = getframe(fig); % Capture t spot.
%         writeVideo(mv,frame); % Save
        errs{j,1}=err; errs{j,2}=ns; errs{j,3}=t(i); errs{j,4}=h1err;
        j=j+1;
    end
%     %%% Write Movie.
%     %
%     close(mv);
end