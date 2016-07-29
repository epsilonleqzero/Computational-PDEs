%% Finite Element Tests
%
%   Written by Ted Kwan for Math 226C Project 1
%
%   This script runs the tests for the finite element method,
%   then plots the convergence results.

clear all;
alpha=@(x,y,t) -25.*(((x-t+0.5).^2)+((y-t+0.5).^2));
alphap=@(x,y,t) ((x-t+0.5)+(y-t+0.5));
gammat=@(t) -100.*((t-0.5).^2);
betat=@(t) (0.1).*(1-exp(gammat(t)));
utrue=@(x,y,t) betat(t).*exp(alpha(x,y,t));
uexact=@(x,t) betat(t).*exp(alpha(x(:,1),x(:,2),t));
u0=@(x,y) utrue(x,y,0);
f=@(x,y,t)  exp(alpha(x,y,t)).*((20.*(t-0.5).*exp(gammat(t)))+ ...
            50.*betat(t).*(2+(x-t+0.5)-50.*((x-t+0.5).^2) + ...
            (y-t+0.5)-50.*((y-t+0.5).^2)));
%% Space/Time Conditions.
%
N=5;
dt=0.2; t0=0; tf=2; h=0.25; % Both space and time
%dt=0.2; t0=0; tf=2; h=0.1; % Space Error
%dt=0.1; t0=0; tf=2; h=0.05; % Time Error
H1errs=zeros(N,1); Ns=zeros(N,1);
L2errs=zeros(N,1); Ns1=zeros(N,1);
[node,elem] = squaremesh([0,1,0,1],h);
for i=1:N
    [ut,node,elem,t,M,A]=FiniteElemHeat(node,elem,f,u0,t0,tf,dt,utrue);
    Nt=length(t); errs=zeros(Nt,1);
    utestL2=utrue(node(:,1),node(:,2),t(end));
    L2errs(i)=sqrt((((utestL2-ut))')*M*(utestL2-ut));
    H1errs(i)=sqrt(((utestL2-ut)')*A*(utestL2-ut));
    Ns(i)=length(node(:,1))+1/(dt);
    [node,elem]=uniformrefine(node,elem);
    dt=dt/2;
end
figure; showrate2(Ns,L2errs,[],'-*','||u_I-u_h||_{0,2}',...
          Ns,H1errs,[],'k-.','|u_I-u_h|_{1,2}');