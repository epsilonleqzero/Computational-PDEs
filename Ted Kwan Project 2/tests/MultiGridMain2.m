%% Multi-Grid Methods
%
%
clear all;
h=1/64;
ftest=@(x,y) (x-x+1); % Silly modification needed to vectorize
%fteste=@(x,y) (1);
%ftrue=@(x,y) fourierPoisson(x,y);
ug=@(x,y) (0);
% %% Finite Difference - Part 1
% %
% %
% %set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.4,0.3]);
%[nodec,elemc] = squarequadmesh([0,1,0,1],h);
%subplot(1,2,1); showmesh(nodec,elemc);
% [nodef,elemf] = uniformrefinequad(nodec,elemc);
% %subplot(1,2,2); showmesh(nodef,elemf);
%[u,x,y,errs,stop]=PoissonFinDifMG(h,ftest,ug,'twogrid',5000,2);
 
% % 
% figure;
% surf(x,y,u); axis tight;
% xlabel('x'); ylabel('y');
% zlabel('u');
% % ut=arrayfun(ftrue,x,y);
% % err=abs(ut-u);

% %% Finite Element - Part 2
% %
% %
% %set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.4,0.3]);

% times=zeros(4,1);
itrs=zeros(4,1);
m=[2:n]; itrs=zeros(n-1,1); hs=zeros(n-1,1);
hs(m-1)=((1/2).^m);
% for J=3:6
     %[node,elem] = circlemesh(0,0,1,h);
%     tic;
     %[u,node,elem,itrs,res]=FiniteElemMG(node,elem,fteste,4,ug,1000,h,'pcg');
        [u,x,y,errs,itrs]=PoissonFinDifMG(h,ftest,ug,'multigrid',10000,6);
%     [u,node,elem,itrs(J-2)]=FiniteElemPCG(node,elem,ftest,J,ug,1000,h);
%     times(J-2)=toc
%     clear u node elem;
% end
% figure;
% plot([3:6],itrs);
% xlabel('J'); ylabel('Iterations Needed');
% title('Iterations J-level V-cycle PCG');
% axis tight;
% grid on;
% figure;
% plot([3:6],times);
% xlabel('J'); ylabel('t - seconds');
% title('Time for J-level V-cycle PCG');
% axis tight;
% grid on;


%showresult(node,elem,u);
% figure; n=[1:length(errs)];
% semilogy(n,errs);
% xlabel('Step'); ylabel('Residual');
% title('Residual Using Gauss-Seidel');
%  axis tight; %grid on;
