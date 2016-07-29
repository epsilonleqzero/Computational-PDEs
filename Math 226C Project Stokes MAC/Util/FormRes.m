%% Form Residual Equation
%
%   Written by Ted Kwan for Math 226C Project 2
%   
%   This function forms the residuals for each component
%   in the Stokes problem. The residuals are formed with respect
%   to the the Marker and Cell (MAC) scheme.
%
function [Resu,Resv,Resp] = FormRes(u,v,p,f1,f2,f3,h)
    %% Initial Setup
    %
    [n,m]=size(p);
    [Au,Av]=AuS(u,v,h);
    Resu=zeros(n,m+1);
    Resv=zeros(n+1,m);
    Resp=zeros(n,m);
    %% U Residual
    %
    i = 1:n; j = 2:m;
    px=(p(i,j)-p(i,j-1))/h;
    Resu(i,j)=-Au(i,j)-px+f1(i,j);
    %% V Residual
    %
    i = 2:n; j = 1:m;
    py=(p(i-1,j)-p(i,j))/h;
    Resv(i,j)=-Av(i,j)-py+f2(i,j);
    %% Continuity Residual
    %
    i=[1:n]; j=[1:m];
    Resp(i,j)=f3(i,j)-(((u(i,j+1)-u(i,j))-(v(i+1,j)-v(i,j)))/h);
end