%% Matrix Free Calculation of Au
%
%   Written by Ted Kwan for Math 226B
%
%   This function is the matrix free version
%   of A*u. Created using code from FDMcode
%   document by Professor Chen.
function [Au] = AuP(u)
    [n,m]=size(u); i=2:n-1; j=2:m-1;
    Au=zeros(n,m);
    Au(i,j) = 4*u(i,j) - u(i-1,j) - u(i+1,j) ...
            - u(i,j-1) - u(i,j+1);
end