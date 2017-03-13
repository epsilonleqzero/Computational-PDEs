%% Poisson Equation Finite Difference
%
%   Written by Ted Kwan For Math 226A
%   Project 2
%   This solves the Poisson equation with Dirichlet
%   boundary conditions.
%
%   The matrix solving method can be chosen by inputting
%   any of the following:
%   * 'Gauss' - Gauss-Seidel.
%   * 'Gaussrb' - Gauss-Seidel Symmetric.
%   * 'Jacobi' - Weighted Jacobi Method.
%   * 'TwoGrid' - Two Grid Method.
%   * 'MultiGrid' - Multi-Grid method (must specify steps)
function [u,x,y,err,stop] = PoissonFinDifMG(h,f,g,type,itrmax ...
                                            ,steps)
    %% Initial Setup
    %
    [x,y]=ndgrid(0:h:1,0:h:1); nx=length(x(1,:)); nxint=nx-2;
    ny=length(y(:,1)); nyint=ny-2;
    dx=h; dy=h; N=(ny)*(nx); u=zeros(N,1);
    isbd = true(ny,nx); isbd(2:end-1,2:end-1) = false;
    isfree=~isbd; bdidx = find(isbd(:));
    intidx=find(isfree(:)); u(bdidx)=g(x(bdidx),y(bdidx));
    %% Setup Residual
    %
    u=reshape(u,ny,nx); fu=(h.^2)*f(x,y);
    %% Initial guess for Gauss-Sidel
    %
    errgr=zeros(size(u)); err=zeros(itrmax,1);
    u(intidx)=ones(size(u(intidx)));
    icur=itrmax; uapx=Psolve(h,fu);
    tol=(h.^2)*max(abs(uapx(:))); 
    errcur=2*tol; stop=1; omega=0.43;
    er=zeros(size(u)); hc=h;
    %% Solve System
    %
    %   Check iterative method used and solves
    %   the linear system.
    if(strcmpi(type,'gauss'))
        %%% Gauss-Seidel
        %
        while (errcur > tol && icur >stop)
            u = GaussSid(u,fu);
            errcur = max(abs(uapx(:)-u(:)));
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    elseif(strcmpi(type,'gaussrb'))
        %%% Gauss-Seidel Symmetric
        %
        while (errcur > tol && icur >stop)
            u = GaussSidRB(u,fu);
            errcur = max(abs(uapx(:)-u(:)));
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    elseif(strcmpi(type,'jacobi'))
        %%% Weighted Jacobi
        %
        while (errcur > tol && icur >stop)
            u=wJac(u,fu,omega);
            errcur = max(abs(uapx(:)-u(:)));
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    elseif(strcmpi(type,'twogrid'))
        %%% Two-Grid
        %
        while (errcur > tol && icur >stop);
            u=TwoGridP(u,fu,h,3,35);
            errcur = max(abs(uapx(:)-u(:)))
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    elseif(strcmpi(type,'multigrid'))
        %%% Multigrid
        %
        while (errcur > tol && icur >stop)
            [u,hc]=MultiGridP(u,fu,h,3,steps);
            errcur = max(abs(uapx(:)-u(:)))
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    else
        %%% Implements mldivide.
        %
        u=Psolve(h,fu); u=reshape(u,ny,nx);
        err=0;
    end
end