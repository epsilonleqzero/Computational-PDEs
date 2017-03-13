%% Classical Iterative Method Tests
%
%

clear all;
soltype={'jacobi','gauss','gaussrb'};
titles={'Jacobi - One Step','Jacobi - Two Step','Jacobi - Three Steps', ...
        'Gauss-Seidel - One Step','Gauss-Seidel - Two Step','Gauss-Seidel - Three Steps'};
ftest=@(x,y) (x-x+1);
ftrue=@(x,y) fourierPoisson(x,y);
ug=@(x,y) (0);
%PoissonFinDifMG(h,f,g,type,tol,us)
% n=7;
% m=[2:n]; itrs=zeros(n-2,1); hs=zeros(n-1,1);
% hs(m-1)=((1/2).^m); times=zeros(n-2,1);
% tols=zeros(n-2,1); hcs=zeros(n-2,1);
% js=[2:n-1]'; smooths=3*ones(n-2,1);
% Ns(m-1)=(1./(hs))+1;
% for i=1:2
%     for j=2:n-1
%     tic
%     [u,x,y,errs,itrs(j-1),tols(j-1),hcs(j-1)]=PoissonFinDifMG(hs(j),ftest,ug,'multigrid',100000,j);
%     times(j-1)=toc;
%     end
%     clear u x y;
% end
% figure; loglog(1./hs(2:end),times); grid on;
% xlabel('1/h'); ylabel('t - Seconds');
% title('MultiGrid Method CPU Time'); axis tight;
% figure; plot(Ns(2:end),itrs); grid on;
% xlabel('Matrix Size N'); ylabel('iterations');
% title('MultiGrid Method Iterations'); axis([8,128,3.8,6.1]);
% mgt=[js,hcs,smooths,tols];
% dlmwrite('mgtests.txt', mgt , '&');
%for i=1:2
    %for j=1:3
     [u,x,y,errgr,itrs(1)]=PoissonFinDifMGErrs((1/16),ftest,ug,soltype{2},3);
%      figure; surf(x,y,errgr{1});
%      title(titles{1});
%      figure; surf(x,y,errgr{2});
%      title(titles{2});
%      figure; surf(x,y,errgr{3});
%      axis([0,1,0,1,0,1])
%     title(titles{3});
%     [u,x,y,errgr,itrs(1)]=PoissonFinDifMGErrs((1/64),ftest,ug,soltype{1},3);
%     figure; 
    subplot(2,2,1); surf(x,y,errgr{1}); title(titles{4});
    subplot(2,2,2); surf(x,y,errgr{2}); title(titles{5});
    subplot(2,2,3); surf(x,y,errgr{3}); title(titles{6}); axis([0,1,0,1,0,1])
    
%end
