%% Multi-Grid Methods
%
%
clear all;
%%% Initial Setup
%
h=1/8;
ftest=@(x,y) (x-x+1); % Silly fixup to vectorize.
ug=@(x,y) (0);
%% Finite Element - Iterations and Time
%
js=[3:6]'; n=length(js);
times=zeros(n,1); itrs=zeros(n,1); tols=zeros(n,1);
for J=3:6
    [node,elem] = circlemesh(0,0,1,h); tic;
    [u,elem,node,itrs(J-2),res,tols(J-2)]=FiniteElemMG(node,elem,ftest ...
                                          ,J,ug,2500,1e-6);
    times(J-2)=toc;
    clear u node elem;
end
%% Save Params
% mgt=[js,itrs,tols,times];
% dlmwrite('mgtests2.txt', mgt , '&');
%%% Plot Iterations and Time
%
figure;
plot([3:6],itrs);
xlabel('J'); ylabel('Iterations Needed');
title('Iterations J-level V-cycle PCG');
axis([3,6,0,20]);
grid on; figure;
plot([3:6],times);
xlabel('J'); ylabel('times');
title('CPU Time Used for PCG');
axis tight; grid on;
% %% Plot Residual
% %
% [node,elem] = circlemesh(0,0,1,h);
% [u,elem,node,~,res,~]=FiniteElemMG(node,elem,ftest ...
%                                            ,4,ug,2500,1e-6);
% figure; n=[1:length(res)];
% plot(n,res);
% xlabel('Steps'); ylabel('Residual');
% title('Residual Using a 4-level V-cycle PCG');
% grid on;
% axis tight;
