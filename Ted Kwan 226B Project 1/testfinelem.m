%% Test Finite Element Method
%
%   Written by Ted Kwan for Math 226B Project 1.
%
%   This checks the error and plots the convergence rate for both
%   the PDE with either Dirichlet or Neumann boundary conditions
%   given in the project.
%
clear all;

%% Dirichlet Boundary Setup
%

% mesht='square';
% bdryt='dirichlet';
% f=@(x,y) (8*(pi^2)*sin(2*pi*x).*cos(2*pi*y));
% uHi=inline('sin(2*pi*pxy(:,1)).*cos(2*pi*pxy(:,2))','pxy');
% Du = inline('[2*pi*cos(2*pi*pxy(:,1)).*cos(2*pi*pxy(:,2)),-2*pi*sin(2*pi*pxy(:,1)).*sin(2*pi*pxy(:,2))]','pxy');
% squ=[0,1,0,1]; h=0.25;
% g=inline('sin(2*pi*x).*cos(2*pi*y)');
% gn=@(x) NeumBdry(x);

%% Neumann Problem Setup
%

mesht='circle';
bdryt='Neumann';
f=@(x,y) (1); g=@(x,y) (-0.5);
DuN=inline('[-pxy(:,1)/2,-pxy(:,2)/2]','pxy');
uHi=inline('(-(pxy(:,1).^2)-(pxy(:,2).^2))/4','pxy');
x=0;y=0;r=1;h=0.2;

%% Refine and Plot
%
%   This section refines the grid and plots the different
%   errors after recording the amount of elements used
%   in the particular refinement.
%
N=5; err=zeros(N,3); ns=zeros(N,1);
for ref=0:N
%    [us,node,elem,A]=FiniteElem(mesht,bdryt,f,ref,g,0,squ,h);
   [us,node,elem,A]=FiniteElem(mesht,bdryt,f,ref,g,g,x,y,r,h);
   ns(ref+1)=length(elem(:,1)); uh1 = uHi(node);
   err(ref+1,3) = max(abs((A*uh1)-(A*us)));  clear A; %Retrieve memory.
   err(ref+1,1) = getH1error(node,elem,DuN,us);
   err(ref+1,2) = getL2error(node,elem,uHi,us);
end
%% Plot Convergence Rate
%
%   This plots the convergence rate with the showrate function
%   given in the ifem folder.
%
sty1='||u-u_h||_H'; sty2='||u-u_h||_{L^{2}}'; sty3='||Au-Au_h||';
r1=showrate(ns,err(:,1),[],'r-+'); r2=showrate(ns,err(:,2),[],'b-+');
r3=showrate(ns,err(:,3),[],'g-+');
%%% Fixup for older versions on Matlab.
%
%   The round command does not work with a second argument
%   in older versions on Matlab, but it can be re-implemented 
%   by just rounding the numbers manually.
%
r1nd=round((100*r1))/100; r2nd=round((100*r2))/100;
r3nd=round((100*r3))/100;
h_legend = legend(sty1,['C_1N^{' num2str(r1nd) '}'],...
                  sty2,['C_2N^{' num2str(r2nd) '}'],...
                  sty3,['C_3N^{' num2str(r3nd) '}'],'LOCATION','Best');
set(h_legend,'FontSize',12);