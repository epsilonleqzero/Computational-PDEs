%% Finite Element Preconditioned Conjugate Gradient
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the preconditioned
%   conjugate gradient method using a J-level
%   v-cycle as the preconditioner. This method is
%   used with the finite element method in order to 
%   solve the system Au=b where a is the stiffness matrix
function [u,k,res] = PCGvCycle(b,u,J,As,Res,Pro,R,freeNodes, ...
                                BdNodes,tol,maxitr)
    k=1; tol=tol; r=b-As{J}*u; r(BdNodes{J+1})=0; 
    res=zeros(maxitr,1); err=2*tol;
    while(err>tol && k<maxitr)
        Br=FEMvCycle(r,J,As,Res,Pro,R,freeNodes,BdNodes);
        rho=r'*Br;
        if(k==1)
            p=Br;
        else
            beta=rho/rho_old; p=Br+beta.*p;
        end
        Ap=As{J}*p; alpha=rho/(p'*As{J}*p);
        u(freeNodes{J+1})=u(freeNodes{J+1})+ ...
                           alpha*p(freeNodes{J+1});
        r=r-alpha*Ap; r(BdNodes{J+1})=0;
        rho_old=rho; k=k+1; err=norm(r)
        res(k)=err;
    end
    res=res(find(res));
end