%% Graph Solutions
%
%   Written by Ted Kwan for Math 226B Project 1.
%
%   This plots the solution to either of the PDEs
%   given in the project description.
%% Setup initial properties
%
%%% Dirichlet Problem on the unit square
%
% mesht='square';
% bdryt='dirichlet';
% f=@(x,y) (8*(pi^2)*sin(2*pi*x).*cos(2*pi*y));
% squ=[0,1,0,1]; h=0.25;
% g=inline('sin(2*pi*x).*cos(2*pi*y)');
%%% Neumann Problem on the unit disk
%
mesht='circle';
bdryt='Neumann';
g=@(x,y) (-0.5);
f=@(x,y) (1);
x=0;y=0;r=1;h=0.2;
%% Find Approximation
%
%[us,node,elem,A]=FiniteElem(mesht,bdryt,f,0,g,gn,squ,h);
[us,node,elem,A]=FiniteElem(mesht,bdryt,f,0,g,g,x,y,r,h);
%% Plot
%
showresult(node,elem,u);