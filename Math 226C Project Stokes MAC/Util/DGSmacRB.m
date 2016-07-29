%% Distributed Gauss-Sidel Relaxation For MAC Scheme
%
%   Written by Ted Kwan for Math 226C
%   
%   This function implements the Distributed Gauss-Seidel
%   method for the Stokes equations to be used alone, or
%   with a multigrid solver.
%
%   This function is specific to the Marker and Cell (MAC)
%   scheme, so it can only run on the staggered grids.
%
function [u,v,p] = DGSmacRB(u,v,p,f1,f2,gp,h,omega)
    [n,m] = size(p);
    %% Update U
    %
    u=GaussSidX(u,h*p,(h^2)*f1,omega);
    %% Update V
    %
    v=GaussSidY(v,h*p,(h^2)*f2,omega);
    %% Calculate Divergence.
    %
    i=[1:n]; j=[1:m];
    rc=((((u(i,j+1)-u(i,j))-(v(i+1,j)-v(i,j)))/h)-gp(i,j));
    rc=rc-mean(rc(:))*ones(size(rc));
    %% DQ Equation
    % Corner points
    i = 2:m-1; j = 2:n-1;
    %% Continuity Equation Residual
    %
    %   This calculates the second iteration of the
    %   Poisson equation for the correction.
    %
    dq=zeros(size(p));
    dq(i,1)=((h^2)*rc(i,1)+dq(i,2)+dq(i+1,1)+dq(i-1,1))/3;
    dq(i,n)=((h^2)*rc(i,n)+dq(i,n-1)+dq(i+1,n)+dq(i-1,n))/3;
    dq(1,j)=((h^2)*rc(1,j)+dq(2,j)+dq(1,j+1)+dq(1,j-1))/3;
    dq(m,j)=((h^2)*rc(m,j)+dq(m-1,j)+dq(m,j+1)+dq(m,j-1))/3;
    dq(1,1)=((h^2)*rc(1,1)+dq(2,1)+dq(1,2))/2;
    dq(m,n)=((h^2)*rc(m,n)+dq(m-1,n)+dq(m,n-1))/2;
    dq(1,n)=((h^2)*rc(1,n)+dq(1,n-1)+dq(2,n))/2;
    dq(m,1)=((h^2)*rc(m,1)+dq(m-1,1)+dq(m,2))/2;
    for i = 2:m-1
        for j = 2:n-1
            dq(i,j)=((h^2)*rc(i,j)+dq(i-1,j)+dq(i+1,j)...
                    +dq(i,j-1)+dq(i,j+1))/4;
        end
    end
    %% Update P
    %
    ip=[2:n-1]; jp=[2:m-1];
    %%% Corner Nodes
    %
    p(1,1)=p(1,1)-(2*dq(1,1)-dq(2,1)-dq(1,2))/(h^2);
    p(m,n)=p(m,n)-(2*dq(m,n)-dq(m-1,n)-dq(m,n-1))/(h^2);
    p(m,1)=p(m,1)-(2*dq(m,1)-dq(m-1,1)-dq(m,2))/(h^2);
    p(1,n)=p(1,n)-(2*dq(1,n)-dq(1,n-1)-dq(2,n))/(h^2);
    %%% Edge Nodes
    %
    p(1,jp)=p(1,jp)-(3*dq(1,jp)-dq(2,jp)-dq(1,jp+1)-dq(1,jp-1))/(h^2);
    p(m,jp)=p(m,jp)-(3*dq(m,jp)-dq(m-1,jp)-dq(m,jp+1)-dq(m,jp-1))/(h^2);
    p(ip,1)=p(ip,1)-(3*dq(ip,1)-dq(ip,2)-dq(ip+1,1)-dq(ip-1,1))/(h^2);
    p(ip,n)=p(ip,n)-(3*dq(ip,n)-dq(ip,n-1)-dq(ip+1,n)-dq(ip-1,n))/(h^2);
    %%% Interior nodes
    %
    p(ip,jp)=p(ip,jp)-(4*dq(ip,jp)-dq(ip-1,jp)-dq(ip+1,jp)...
            -dq(ip,jp-1)-dq(ip,jp+1))/(h^2);
    p=p-mean(p(:))*ones(size(p));

    %% Apply Correction to U
    %
    i=[1:n]; j=[2:m];
    u(i,j)=u(i,j)+(dq(i,j)-dq(i,j-1))/h;
    %% Apply Correction to V
    %
    i=[2:n]; j=[1:m];
    v(i,j)=v(i,j)+(dq(i-1,j)-dq(i,j))/h;
end