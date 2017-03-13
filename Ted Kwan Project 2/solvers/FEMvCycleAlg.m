%% Algebraic Multi-Grid V-Cycle
%
%   Written by Ted Kwan for Math 226B
%
%   This function implements a J-level v-cycle to be used alone
%   or as the preconditioner for the conjugate gradient method.
%   This function is used with the finite element method in
%   order to solve the system Au=b where A is the stiffness
%   matrix
function [e] = FEMvCycleAlg(r,J,As,Res,Pro,R,mu)
    ri=cell(J,1); ei=cell(J,1); ri{J}=r;
    for i=J:-1:2
        %%% Pre Smoothing
        %
        ei{i}=R{i,1}\(ri{i});
         for j=1:(mu-1) % Fixup because first step
                        % is outside of loop
             ei{i}=ei{i}+R{i,1}\(ri{i}-As{i}*ei{i});
         end
        %%% Update and Restrict
        %
        ri{i-1}=Res{i}*(ri{i}-As{i}*ei{i});
    end
    %%% Solve on Coarse Grid
    %
    ei{1}=As{1}\ri{1};
    for i=2:J
        %%% Prolongate and Correct
        %
        ei{i}=ei{i}+Pro{i}*ei{i-1};
        %%% Post Smoothing
        %
        for j=1:mu
            ei{i}=ei{i}+R{i,2}\(ri{i}-As{i}*ei{i});
        end
    end
    e=ei{J};
end