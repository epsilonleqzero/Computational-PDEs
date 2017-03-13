%% Create Matrix and solve
%
%   Written by Ted Kwan for Math 226B
%
%   This method solves the matrix equation for
%   the finite difference method.
function [u] = Psolve(h,fu)
    %%% Setup System
    %
    b=fu(:); nx=(1/h)+1; ny=(1/h)+1;
    N=(ny)*(nx); u=zeros(N,1);
    isbd = true(ny,nx);
    isbd(2:end-1,2:end-1) = false;
    isfree=~isbd; intidx=find(isfree(:));
    %%% Create Matrix
    %
    ex = ones(nx,1); ey = ones(ny,1);
    Tx = spdiags([-ex 2*ex -ex], -1:1, nx, nx);
    Ty = spdiags([-ey 2*ey -ey], -1:1, ny, ny);
    A = kron(speye(nx),Ty) + kron(Tx,speye(ny));
    %%% Solve
    %
    u(intidx)=A(intidx,intidx)\b(intidx); u=reshape(u,ny,nx);
end