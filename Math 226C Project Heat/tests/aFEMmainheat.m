%% Adaptive Finite Element Tests
%
%   Written by Ted Kwan for Math 226C Project 1
%
%   This script runs the tests for the adaptive finite
%   element method, then plots the convergence results.

clear all;
%% Initial Setup
%
alpha=@(x,y,t) -25.*(((x-t+0.5).^2)+((y-t+0.5).^2));
alphap=@(x,y,t) ((x-t+0.5)+(y-t+0.5));
gammat=@(t) -100.*((t-0.5).^2);
betat=@(t) (0.1).*(1-exp(gammat(t)));
utrue=@(x,y,t) betat(t).*exp(alpha(x,y,t));
uexact=@(x,t) betat(t).*exp(alpha(x(:,1),x(:,2),t));

f=@(x,y,t)  exp(alpha(x,y,t)).*((20.*(t-0.5).*exp(gammat(t)))+ ...
            50.*betat(t).*(2+(x-t+0.5)-50.*((x-t+0.5).^2) + ...
            (y-t+0.5)-50.*((y-t+0.5).^2)));
%dt=0.1; t0=0; tf=1.4; h=0.05; %Use this for H1 error.
%dt=0.025; t0=0; tf=1.4; h=0.05; %Use this for H1 error.
%dt=0.025; t0=0; tf=1.8; h=0.05; %Use this for L2 error.
dt=0.1; t0=0; tf=1.8; h=0.2; %Use this for L2 error.
u0=@(x,y) utrue(x,y,0);
h2=0.01;
Du=@(x) -50.*betat(tf).*[exp(alpha(x(:,1),x(:,2),tf)).*(x(:,1)-tf+0.5), ...
          exp(alpha(x(:,1),x(:,2),tf)).*(x(:,2)-tf+0.5)];
laplaceu=@(x) 50.*betat(tf).*exp(alpha(x(:,1),x(:,2),tf)).*(-2+...
          50.*((x(:,1)-tf+0.5).^2)+50.*((x(:,2)-tf+0.5).^2));
[nodestart,elemstart] = squaremesh([0,1,0,1],h);
[nodestart1,elemstart1] = squaremesh([0,1,0,1],h2);

%% Finite Element Method
%
[ut1,node1,elem1,t1,M1,A2,At2,L2errs]=afemheat(nodestart,elemstart,f,u0,t0,tf,dt,utrue,1,1);
%[ut1,node1,elem1,t1,M1,A2,At2,L2errs]=astfemheat(nodestart,elemstart,f,u0,t0,tf,dt,utrue,1);
[ut4,node5,elem5,t5,M5,A5,At5]=FiniteElemHeat(nodestart1,elemstart1,f,u0,t0,tf,dt,utrue);

%% Record Errors and plot
%
N1=L2errs{1,2};
aL2err=L2errs{1,1};
aH1err=L2errs{1,4};
figure;
showrate2(N1,aL2err,[],'-*','||u_I-u_h||_{0,2}', ...
          N1,aH1err,[],'-*','||u_I-u_h||_{1,2}');
N2=L2errs{2,2};
aL2err2=L2errs{2,1};
aH1err2=L2errs{2,4};
figure;
showrate2(N2,aL2err2,[],'-*','||u_I-u_h||_{0,2}', ...
          N2,aH1err2,[],'-*','||u_I-u_h||_{1,2}');
N3=L2errs{3,2};
aL2err3=L2errs{3,1};
aH1err3=L2errs{3,4};
figure;
showrate2(N3,aL2err3,[],'-*','||u_I-u_h||_{0,2}', ...
          N3,aH1err3,[],'-*','||u_I-u_h||_{1,2}');

utestL21=utrue(node1(:,1),node1(:,2),t1(end));
L2err1=sqrt((((utestL21-ut1))')*M1*(utestL21-ut1))
H1err1=sqrt((((utestL21-ut1))')*A2*(utestL21-ut1))

utestL24=utrue(node5(:,1),node5(:,2),t5(end));
L2err4=sqrt((((utestL24-ut4))')*M5*(utestL24-ut4))
H1err4=sqrt(((utestL24-ut4)')*A5*(utestL24-ut4))
