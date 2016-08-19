%% Multi-Grid Methods
%
%
clear all;
%%% Initial Setup
%
h=1/8;
ub=@(x,y) (0.1+((x+2*y)./sqrt(5)));
utrue=@(x,y) -log((1+cos(ub(x,y)))./(1-cos(ub(x,y)))); %Silly fixup to vectorize.
f=@(x,y) ((utrue(x,y)));
ug=@(x,y) (0);
%% Finite Element - Iterations and Time
%
js=[3:6]'; n=length(js);
%times=zeros(n,1); itrs=zeros(n,1); tols=zeros(n,1);
%for J=3:6
    [node,elem] = squaremesh([0,1,0,1],h); %tic;
    [u,node,elem]=FiniteElemMG(node,elem,f ...
                                         ,4,utrue,150,1e-9);
    %times(J-2)=toc;
    showsolution(node,elem,u);
    %clear u node elem;
%end
%[x,y]=meshgrid()
[x,y]=meshgrid(0:0.01:1,1:-0.01:0);
usol=utrue(x,y);
figure;
surf(x,y,usol);
%% Save Params
% mgt=[js,itrs,tols,times];
% dlmwrite('mgtests2.txt', mgt , '&');
%%% Plot Iterations and Time
%
% figure;
% plot([3:6],itrs);
% xlabel('J'); ylabel('Iterations Needed');
% title('Iterations J-level V-cycle PCG');
% axis([3,6,0,20]);
% grid on; figure;
% plot([3:6],times);
% xlabel('J'); ylabel('times');
% title('CPU Time Used for PCG');
% axis tight; grid on;
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
