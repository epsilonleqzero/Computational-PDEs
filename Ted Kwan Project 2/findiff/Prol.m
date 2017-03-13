%% Prolongation
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the bilinear
%   prolongation map for the multi-grid
%   method.
function [uf] = Prol(uc,Nc)
    %%% Initial Setup
    %
    Nf=2*Nc-1; jc=[2:Nc-1];
    jf=[2.*jc-1]; uf=zeros(Nf,Nf);
    %%% Center
    %
    uf(jf,jf)=uc(jc,jc);
    %%% Sides
    %
    uf(jf-1,jf)=uf(jf-1,jf)+(0.5*uc(jc,jc));
    uf(jf,jf-1)=uf(jf,jf-1)+(0.5*uc(jc,jc));
    uf(jf+1,jf)=uf(jf+1,jf)+(0.5*uc(jc,jc));
    uf(jf,jf+1)=uf(jf,jf+1)+(0.5*uc(jc,jc));
    %%% Corners
    %
    uf(jf+1,jf+1)=(uf(jf+1,jf+1)+(0.25*uc(jc,jc)));
    uf(jf+1,jf-1)=(uf(jf+1,jf-1)+(0.25*uc(jc,jc)));
    uf(jf-1,jf+1)=(uf(jf-1,jf+1)+(0.25*uc(jc,jc)));
    uf(jf-1,jf-1)=(uf(jf-1,jf-1)+(0.25*uc(jc,jc)));
end