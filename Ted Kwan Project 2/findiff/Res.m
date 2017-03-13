%% Restriction
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the bilinear
%   restriction map for the multi-grid
%   method.
function [uc] = Res(uf,Nc)
    Nf=2*Nc-1; jc=[2:Nc-1];
    jf=[2.*jc-1]; uc=zeros(Nc,Nc);
    uc(jc,jc)=(4*uf(jf,jf)+2*(uf(jf-1,jf)+uf(jf,jf-1)+...
            uf(jf+1,jf)+uf(jf,jf+1))+uf(jf-1,jf+1)+...
            uf(jf+1,jf-1)+uf(jf-1,jf-1)+uf(jf+1,jf+1))./4;
end