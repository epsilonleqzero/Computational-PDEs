%% Matrix Free Calculation of Au
%
%   Written by Ted Kwan for Math 226C
%
%   This function is the matrix free version
%   of Au and Av for the Stokes Equation using the 
%   Marker and Cell (MAC) scheme. The boundary conditions
%   are applied within the function, so they are not
%   included in the calculation.
%
%   Created using code from FDMcode document by 
%   Professor Chen.
function [Au,Av] = AuS(u,v,h)
    %% Calculate Au
    %
    [nu,mu]=size(u); i=2:nu-1; j=2:mu-1;
    Au=zeros(nu,mu);
    Au(i,j)=(4*u(i,j)-u(i-1,j)-u(i+1,j) ...
            -u(i,j-1)-u(i,j+1))/(h^2);
    Au(nu,j)=(6*u(nu,j)-(4/3)*u(nu-1,j)-u(nu,j-1)-u(nu,j+1))/(h^2);
    Au(1,j)=(6*u(1,j)-(4/3)*u(2,j)-u(1,j-1)-u(1,j+1))/(h^2);
    %% Calculate Av
    %
    [nv,mv]=size(v); i=2:nv-1; j=2:mv-1;
    Av=zeros(nv,mv);
    Av(i,j)=(4*v(i,j)-v(i-1,j)-v(i+1,j) ...
            -v(i,j-1)-v(i,j+1))/(h^2);
    Av(i,mv)=(6*v(i,mv)-(4/3)*v(i,mv-1)-v(i-1,mv)-v(i+1,mv))/(h^2);
    Av(i,1)=(6*v(i,1)-(4/3)*v(i,2)-v(i-1,1)-v(i+1,1))/(h^2);
end