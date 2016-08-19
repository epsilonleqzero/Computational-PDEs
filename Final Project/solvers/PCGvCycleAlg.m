%% Finite Element Preconditioned Conjugate Gradient
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the preconditioned
%   conjugate gradient method using a J-level
%   v-cycle as the preconditioner. This method is
%   used with the finite element method in order to 
%   solve the system Au=b where A is the stiffness matrix
function [u,k,res] = PCGvCycleAlg(b,u,J,As,Res,Pro,R,tol,maxitr,mu)
    %%% Initial Setup
    %
    nmb=norm(b);tol=tol*nmb; k=1; err=2*tol;
    r=b-As{J}*u; res=zeros(maxitr,1);
    %% Main Loop
    %
    while(err>tol && k<maxitr)
        Br=FEMvCycleAlg(r,J,As,Res,Pro,R,mu);
        rho=r'*Br;
        if(k==1)
            p=Br;
        else
            beta=rho/rho_old; p=Br+beta.*p;
        end
        Ap=As{J}*p; alpha=rho/(p'*Ap);
        u=u+alpha*p; r=r-alpha*Ap; rho_old=rho;
        k=k+1; err=norm(r);
        res(k)=err;
    end
    res=res(find(res));
end