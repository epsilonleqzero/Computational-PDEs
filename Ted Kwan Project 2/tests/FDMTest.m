%% Classical Iterative Method Tests
%
%

clear all;
soltype={'jacobi','gauss','gaussrb','twogrid','multigrid'};
titles={'Weighted Jacobi','Gauss-Seidel','Gauss-Seidel Symmetric', ...
            'Two Grid Method','Multi-Grid Method'};

ftest=@(x,y) (x-x+1);
ftrue=@(x,y) fourierPoisson(x,y);
ug=@(x,y) (0);
%PoissonFinDifMG(h,f,g,type,tol,us)
n=7;  hs=zeros(n-1,1);
m=[2:n]; 
hs(m-1)=((1/2).^m);
errs=cell(n,1);
itrs=zeros(n-1,2);
N(m-1)=(1./(hs(m-1)))+1;
for i=1:3
    step=i
    if(i>=4)
        itrs=zeros(n-2,1);
        for j=2:n-1
            [u,x,y,errs{j-1},itrs(j-1)]=PoissonFinDifMG(hs(j) ...
                                ,ftest,ug,soltype{i},100000,j);
        end
        figure;
        plot(1./hs(2:end),itrs); grid on; axis([8,130,0,10]);
        title(titles{i});
    else
        itrs=zeros(n-1,1);
        for j=1:n-1
            [u,x,y,errs{j},itrs(j,i)]=PoissonFinDifMG(hs(j) ...
                                ,ftest,ug,soltype{i},100000,j);
        end
        figure; loglog(1./hs,itrs,'b');
        p = polyfit(log(hs)./log(2),log(itrs(:,i))./log(2),1);
        r = single(p(1));
        rrnd=(round(100*r)/100);
        title(['Rate of convergence is Ch^{' num2str(rrnd) '}'],'FontSize', 14);
        clear r; axis tight;
    end
    xlabel('1/h'); ylabel('iterations'); 
end