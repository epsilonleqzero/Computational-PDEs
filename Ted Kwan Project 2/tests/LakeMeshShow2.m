%% Algebraic Multi-Grid Method Testing
%
%   Written by Ted Kwan for Math 226b
%
%   This script tests the algebraic multi-grid 
%   solver.

%% Initial setup
%
clear all;
load lakemesh.mat;
f=@(x,y) (1); ug=@(x,y) (0);
itrs=zeros(4,3); times=zeros(4,3);
%% Main Loop
%
i=1;
for i=1:3
    N=length(node(:,1)); r=rand(N,1);
    for J=3:6
        tic;
        %%% Main Solve
        %
        % change the last two entries from ug,f to r 
        % for the algebriac PCG method.
        [u,node,elem,itrs(J-2,i)]=AlgebraicPCG(node,elem, ...
                                    J,1000,1e-6,r);
        times(J-2,i)=toc;
    end
    if(i<3)
       [node,elem] = uniformrefine(node,elem);
    end
end
%% Plot Results
%
%%% Iterations
%
figure;
plot([3:6],itrs(:,1),[3:6],itrs(:,2),[3:6],itrs(:,3));
xlabel('J'); ylabel('Iterations Needed');
title('Iterations J-level V-cycle PCG');
axis([3,6,0,40]); grid on;
legend('no refinements','1 refinement','2 refinements' ...
        ,'location','best');
%%% Times
%
figure;
plot([3:6],times(:,1),[3:6],times(:,2),[3:6],times(:,3));
xlabel('J'); ylabel('t - seconds');
title('Time for J-level V-cycle PCG');
axis tight; grid on;
legend('no refinements','1 refinement','2 refinements' ...
        ,'location','best');
%%% Result on Lake Superior
%
%figure; showsolution(node,elem,u);