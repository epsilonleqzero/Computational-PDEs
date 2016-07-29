%% Stokes MAC Scheme
%
%   Written by Ted Kwan For Math 226C Project 3
%   This solves the Stokes Equations with Dirichlet Boundary
%   Conditions using the Marker and Cell (MAC) scheme.
%
%   The matrix solving method can be chosen by inputting
%   any of the following:
%   * 'Gauss' - Gauss-Seidel.
%   * 'TwoGrid' - Two Grid Method.
%   * 'MultiGrid' - Multi-Grid method (must specify steps)
%
function [u,v,p,xp,yp,xu,yu,xv,yv,x,y,err,stop,fu,tol] = StokesMAC(h,f1,f2,g1,g2,g3,itrmax,type,level)
    %% Initial Setup
    %
    [x,y]=meshgrid(-1:h:1,1:-h:-1);
    %%% Pressure Mesh
    %
    [xp,yp]=meshgrid(-1+(h/2):h:1-(h/2),1-(h/2):-h:-1+(h/2));
    %%% U Mesh
    %
    [xu,yu]=meshgrid(-1:h:1,1-(h/2):-h:-1+(h/2));
    %%% V Mesh
    %
    [xv,yv]=meshgrid(-1+(h/2):h:1-(h/2),1:-h:-1);
    
    nyp=length(yp(:,1)); nxp=length(xp(1,:));
    nyu=length(yu(:,1)); nxu=length(xu(1,:));
    nyv=length(yv(:,1)); nxv=length(xv(1,:));
    Nu=(nyu)*(nxu); Nv=(nyv)*(nxv); Np=(nyp)*(nxp);
    %%% Initialize Vectors
    %
    v=zeros(Nv,1); u=zeros(Nu,1);
    p=zeros(Np,1); gp=zeros(Np,1);
    gbdv=zeros(Nv,1); gbdu=zeros(Nu,1);
    %% Boundary and Free-Node Setup
    %
    %%% U
    %
    isbdu = true(nyu,nxu); isbdug = true(nyu,nxu);
    isbdu(1:end,2:end-1) = false; isfreeu=~isbdu;
    bdidxu = find(isbdu(:)); intidxu=find(isfreeu(:));
    %%% V
    %
    isbdv = true(nyv,nxv); isbdv(2:end-1,1:end) = false;
    isfreev=~isbdv; bdidxv = find(isbdv(:));
    intidxv=find(isfreev(:));
    %%% P
    %
    isbdp = true(nyp,nxp); isbdp(2:end-1,2:end-1) = false;
    isfreep=~isbdp; intidxp=find(isfreep(:));
    %% Initialize Interior and BC
    % 
    u(intidxu)=zeros(size(u(intidxu)));
    v(intidxv)=zeros(size(v(intidxv)));
    p(intidxp)=zeros(size(p(intidxp)));
    u(bdidxu)=g1(xu(bdidxu),yu(bdidxu));
    v(bdidxv)=g2(xv(bdidxv),yv(bdidxv));
    gu=g1(x,y); gv=g2(x,y);
    u=reshape(u,nyu,nxu); v=reshape(v,nyv,nxv);
    p=reshape(p,nyp,nxp); gp=reshape(gp,nyp,nxp);
    %% Boundary Conditions for Staggered Boundary
    %
    gbdu=reshape(gbdu,nyu,nxu); gbdv=reshape(gbdv,nyv,nxv);
    j=[2:nxp];
    gbdu(nxp,j)=(8/3)*gu(nxp+1,j)/(h^2);
    gbdu(1,j)=(8/3)*gu(1,j)/(h^2);
    gbdv(j,nxp)=(8/3)*gv(j,nxp+1)/(h^2);
    gbdv(j,1)=(8/3)*gv(j,1)/(h^2);
    fu=f1(xu,yu)+(gbdu); fv=f2(xv,yv)+(gbdv);
    %% Initial guess for Gauss-Sidel
    %
    err=zeros(itrmax,1); icur=itrmax;
    %tol=1e-6;
    [ru,rv,rp]=FormRes(u,v,p,fu,fv,gp,h);
    Ress=[norm(ru),norm(rv),norm(rp)];
    tol=(1e-8)*max(Ress);
    errcur=2*tol; stop=1; 
    %% Solve System
    %
    %   Check iterative method used and solves
    %   the linear system.
    %
    if(strcmpi(type,'mg') || strcmpi(type,'multigrid'))
        %[ru,rv,rp]=FormRes(u,v,p,fu,fv,gp,h);
        %% Run Multi-Grid Method.
        %
        while ((errcur > tol)&& icur >stop)
            %eu=zeros(size(u)); ev=zeros(size(v)); ep=zeros(size(p));
            [u,v,p] = MultiGridStokes(u,v,p,fu,fv,gp,h,4,level);
            [ru,rv,rp]=FormRes(u,v,p,fu,fv,gp,h);
%             [eu,ev,ep] = MultiGridStokes(eu,ev,ep,ru,rv,rp,h,4,level);
%             u=u+eu; v=v+ev; p=p+ep;
            %[ru,rv,rp]=FormRes(u,v,p,fu,fv,gp,h);
            Ress=[norm(ru),norm(rv),norm(rp)];
            errcur = max(Ress)
            err(stop)=errcur; stop=stop+1;
        end

    elseif(strcmpi(type,'twogrid') || strcmpi(type,'2grid'))
        %% Run Two-Grid Method.
        %
        while ((errcur > tol)&& icur >stop)
            [u,v,p] = MultiGridStokes(u,v,p,fu,fv,gp,h,3,2);
            [ru,rv,rp]=FormRes(u,v,p,fu,fv,gp,h);
            Ress=[norm(ru),norm(rv),norm(rp)];
            errcur = max(Ress)
            err(stop)=errcur; stop=stop+1;
        end
    else
        %% Solve with Gauss-Seidel
        %
        omega=2/(1+sin(pi*h));
        while ((errcur > tol )&& icur >stop)
            [u,v,p]=DGSmacRB(u,v,p,fu,fv,gp,h,omega);
            [ru,rv,rp]=FormRes(u,v,p,fu,fv,gp,h);
            errcur = max([norm(ru),norm(rv),norm(rp)])
            err(stop)=errcur; stop=stop+1;
        end
    end
    err=err(find(err));
end