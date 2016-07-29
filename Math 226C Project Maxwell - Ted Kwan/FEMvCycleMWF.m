%% Multi-Grid V-Cycle for Maxwell Equations
%
%   Written by Ted Kwan for Math 226C Project Maxwell
%
%   This function implements a J-level v-cycle 
%   This function is used with the finite element method in
%   order to solve the saddle-point system:
%   
%   [A+B'B, B][u]=[f+B'g]
%   [B',    0][p]=[g]
%
function [eu,ep] = FEMvCycleMWF(ru,rp,J,As,Bs,Ms,Res,Pro,Prop,Resp,R,mu)
    rui=cell(J,1); eui=cell(J,1); rui{J}=ru;
    rpi=cell(J,1); epi=cell(J,1); rpi{J}=rp;
    for i=J:-1:2
        %%% Pre Smoothing
        %
        epi{i}=zeros(size(rpi{i})); eui{i}=zeros(size(rui{i}));
        for j=1:(mu) 
            [eui{i},epi{i}] = DGSMaxwell(rui{i},rpi{i},epi{i},Bs{i},R{i,1},R{i,3});
            gp=-Bs{i}*(eui{i});
        end
        %%% Update and Restrict
        %
        %gp=-Bs{i}*(eui{i});
        rui{i-1}=Res{i}*(rui{i}-As{i}*eui{i}-Bs{i}'*epi{i}+(Ms{i}*Bs{i})'*gp);
        rpi{i-1}=Resp{i}*(rpi{i}+gp);
    end
    %%% Solve on Coarse Grid
    %
    NEint=size(rui{1},1); Nint=size(rpi{1},1);
    Asys = [As{1},Bs{1}';...
            Bs{1},sparse(Nint,Nint)];
    tmpu = Asys\[rui{1}; rpi{1}];
    eui{1}=tmpu(1:NEint);
    epi{1}=tmpu(NEint+1:NEint+Nint);
    for i=2:J
        %%% Prolongate and Correct
        %
        eui{i}=eui{i}+Pro{i}*eui{i-1};
        epi{i}=epi{i}+Prop{i-1}*epi{i-1};
        %%% Post Smoothing
        %
        for j=1:mu
            [eui{i},epi{i}] = DGSMaxwell(rui{i},rpi{i},epi{i},Bs{i},R{i,2},R{i,4});
        end
    end
    eu=eui{J};
    ep=epi{J};
end