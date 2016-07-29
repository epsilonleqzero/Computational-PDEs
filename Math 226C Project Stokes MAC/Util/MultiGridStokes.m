%% Recursive Multigrid V-Cycle
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the recursive
%   multi-grid method for the Stokes Equation
%   equation on the domain $\Omega=[-1,1]\times [-1,1]$
%   with a uniform grid using the Marker and Cell (MAC)
%   scheme.
%
function [u,v,p] = MultiGridStokes(u,v,p,ru,rv,rp,h,m,itr)

    omega=2/(1+sin(pi*h));
    if(itr==1)
    %%% Coarse Grid
    %   Solve Using Gauss-Seidel.
    
    %%% Set Optimal Omega
    %
        icur=5000; tol=1e-8;
        errcurin=2*tol; stop=1; 
        while ((errcurin > tol)&& icur >stop)
            [u,v,p]=DGSmacRB(u,v,p,ru,rv,rp,h,omega);
            [ruf,rvf,rpf]=FormRes(u,v,p,ru,rv,rp,h);
            errcurin = max([norm(ruf),norm(rvf),norm(rpf)]);
            stop=stop+1;
        end
        return;
    end
    %%% Setup index map
    %
    Nc=round(length(p(:,1))/2);
    %%% Pre-Smoothing
    %
    for i=1:m
        [u,v,p]=DGSmacRB(u,v,p,ru,rv,rp,h,1);
    end
    %%% Restriction
    %
    [ruf,rvf,rpf]=FormRes(u,v,p,ru,rv,rp,h);
    [ruc,rvc,rpc]=ResS(ruf,rvf,rpf);
    %%% Coarse Grid Correction
    %
    uc=zeros(Nc,Nc+1); vc=zeros(Nc+1,Nc); pc=zeros(Nc,Nc);
    [eu,ev,ep] = MultiGridStokes(uc,vc,pc,ruc,rvc,rpc,2*h,m,itr-1);
    %%% Post Smoothing
    %
    [eu,ev,ep]=ProS(eu,ev,ep);
    u=u+eu;
    v=v+ev;
    p=p+ep;
    for i=1:m
        [u,v,p]=DGSmacRB(u,v,p,ru,rv,rp,h,1);
    end
end