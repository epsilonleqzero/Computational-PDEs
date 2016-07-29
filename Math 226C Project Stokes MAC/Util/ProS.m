%% Prolongation
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements the bilinear
%   prolongation map for the multi-grid
%   method.
function [uf,vf,pf] = ProS(uc,vc,pc)
    %%% Initial Setup
    %
    Nc=length(pc(:,1)); Nf=2*Nc;
    %uf=zeros(Nfu,Nfu);
    pf=zeros(Nf,Nf);
    %ic=[1:Nc]; jc=[1:Nc];
    %i=[2*jc-1:2*jc]; j=[2.*jc-1];
    %pf(i,j)=pc(ic,jc);
    for i = 1:Nc
        for j = 1:Nc
            pf(2*i-1:2*i,2*j-1:2*j) = pc(i,j);
        end
    end
    %% U
    %
    uf=zeros(Nf,Nf+1); vf=zeros(Nf+1,Nf);
    %%% Red Points
    %
    ic=[2:Nc-1]; jc=[2:Nc];
    i=[2.*ic-1]; j=[2.*jc-1];
    uf(i-1,j)=(3*uc(ic,jc)+uc(ic-1,jc))/4;
    uf(i,j)=(3*uc(ic,jc)+uc(ic+1,jc))/4;
    vf(j,i-1)=(3*vc(jc,ic)+vc(jc,ic-1))/4;
    vf(j,i)=(3*vc(jc,ic)+vc(jc,ic+1))/4;
    %%% Edge Correction
    %
    uf(1,j)=uc(1,jc)/2; uf(Nf,j)=uc(Nc,jc)/2;
    uf(2,j)=(3*uc(2,jc)+uc(3,jc))/4;
    uf(Nf-1,j)=(3*uc(Nc-1,jc)+uc(Nc,jc))/4;
    vf(j,1)=vc(jc,1)/2; vf(j,Nf)=vc(jc,Nc)/2;
    vf(j,2)=(3*vc(jc,2)+vc(jc,2))/4;
    vf(j,Nf-1)=(3*vc(jc,Nc-1)+vc(jc,Nc))/4;
    %%% Black Points
    %
    i=[1:2:Nf];j=[2:2:Nf];
    uf(i,j)=(uf(i,j-1)+uf(i,j+1))/2;
    vf(j,i)=(vf(j-1,i)+vf(j+1,i))/2;
end