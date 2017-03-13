%% Poisson Equation Finite Difference
%
%   Written by Ted Kwan For Math 226A
%   Project 2
%   This solves the Poisson equation with Dirichlet
%   boundary conditions.
%
%   The matrix solving method can be chosen by either
%   inputting 'Gauss' as the type, or any other string.

function [u,x,y,errgr,stop] = PoissonFinDifMGErrs(h,f,g,type,itrmax,steps)
    %%% Initial Setup
    %
    [x,y]=ndgrid(0:h:1,0:h:1);
    nx=length(x(1,:)); nxint=nx-2;
    ny=length(y(:,1)); nyint=ny-2;
    dx=h; dy=h; N=(ny)*(nx); u=zeros(N,1);
    isbd = true(ny,nx); isbd(2:end-1,2:end-1) = false;
    isfree=~isbd; bdidx = find(isbd(:));
    intidx=find(isfree(:)); u(bdidx)=g(x(bdidx),y(bdidx));
    %%% Setup Residual
    %
    u=reshape(u,ny,nx); fu=(h.^2)*f(x,y);
    %r=fu;
    %r=zeros(size(fu));
    %r(intidx)=fu(intidx);
    %Au=AuP(u(intidx)); %res=(b-Au);
    %%% Initial guess for Gauss-Sidel
    %
    errgr=zeros(size(u)); err=zeros(itrmax,1);
    u(intidx)=rand(size(u(intidx)));
    icur=itrmax; uapx=Psolve(h,fu);
    tol=(h.^2)*max(abs(uapx(:))); 
    errcur=2*tol; stop=1; omega=0.43;
    er=zeros(size(u));
    errgr=cell(3,1);
    %%% Checks for Matrix solving type.
    %
    if(strcmpi(type,'gauss'))
        while (errcur > tol && icur >=stop)
            u = GaussSid(u,fu);
            errs=abs(uapx-u);
            errgr{stop}=errs;
            errcur = max(errs(:));
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    elseif(strcmpi(type,'gaussrb'))
        while (errcur > tol && icur >=stop)
            u = GaussSidRB(u,fu);
            errs=abs(uapx-u);
            errgr{stop}=errs;
            errcur = max(errs(:));
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    elseif(strcmpi(type,'jacobi'))
        while (errcur > tol && icur >=stop)
            u=wJac(u,fu,omega);
            errs=abs(uapx-u);
            errgr{stop}=errs;
            errcur = max(errs(:));
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    elseif(strcmpi(type,'twogrid'))
        while (errcur > tol && icur >stop);
            u=TwoGridP(u,fu,h,3,100);
            errcur = max(abs(uapx(:)-u(:)));
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    elseif(strcmpi(type,'multigrid'))
        while (errcur > tol && icur >stop)
            u=MultiGridP(u,fu,h,3,steps);
            errcur = max(abs(uapx(:)-u(:)));
            err(stop)=errcur; stop=stop+1;
        end
        err=err(find(err));
    else
        %%% Implements mldivide.
        %
        fu=(h^2)*f(x(2:end-1,2:end-1),y(2:end-1,2:end-1));
        ex = ones(nxint,1);
        Tx = spdiags([-ex 2*ex -ex], -1:1, nxint, nxint);
        ey = ones(nyint,1);
        Ty = spdiags([-ey 2*ey -ey], -1:1, nyint, nyint);
        A = kron(speye(nxint),Ty) + kron(Tx,speye(nyint));
        u(intidx)=A\fu(:);
        u=reshape(u,ny,nx);
        err=0;
    end
end