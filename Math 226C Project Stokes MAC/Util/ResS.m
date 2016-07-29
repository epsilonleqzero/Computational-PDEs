%% Restriction
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the bilinear
%   restriction map for the multi-grid
%   method.
function [uc,vc,pc] = ResS(uf,vf,pf)
    Nc=round(length(pf(:,1))/2);
    %% Restrict Pressure
    %
    pc=zeros(Nc,Nc);
    jc=[1:Nc]; ic=[1:Nc];
    j=[2.*jc-1]; i=[2.*ic-1];
    pc(ic,jc)=(pf(i,j)+pf(i,j+1)+pf(i+1,j)+pf(i+1,j+1))/4;
    %% Restrict u
    %
    uc=zeros(Nc,Nc+1);
    jc=[2:Nc]; ic=[1:Nc-1];
    j=[2.*jc-1]; i=[2.*ic-1];
    uc(ic,jc) = (uf(i,j-1)+2*uf(i,j)+uf(i,j+1) ...
              +uf(i+1,j-1)+2*uf(i+1,j)+uf(i+1,j+1))/8;
    %% Restrict v
    %
    vc=zeros(Nc+1,Nc);
    jc=[1:Nc-1]; ic=[2:Nc];
    j=[2.*jc-1]; i=[2.*ic-1];
    vc(ic,jc) = (vf(i-1,j)+vf(i-1,j+1) ...
                +2*vf(i,j)+2*vf(i,j+1)...
                +vf(i+1,j)+vf(i+1,j+1))/8;
end