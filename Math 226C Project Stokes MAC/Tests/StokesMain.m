%% Stokes Equation MAC Main File
%
%   Written by Ted Kwan for Math 226C Project 2
%   This script runs the different methods.
%
clear all;
f1=@(x,y) (x-x);
f2=@(x,y) (x-x);
% %% First Problem
% %
% utrue=@(x,y) (20.*(x).*(y.^3));
% vtrue=@(x,y) 5*((x.^4)-(y.^4));
% ptrue=@(x,y) (60*((x.^2).*y)-20*(y.^3));

%% Lid Driven Cavity Problem
%
utrue=@(x,y) (arrayfun(@(x,y) LidCavityu(x,y),x,y));
vtrue=@(x,y) (x-x);
ptrue=@(x,y) (x-x);

%%% V-Cycle Tests
%
%r=[4:6]'; hs=(2^(-7))*ones(length(r),1);
r=[4:7]'; hs=(2.^(-r-1));
%% Two Grid Test
%
%r=[3:6]'; hs=(2.^(-r));

Nh=length(hs); Ns=zeros(length(hs),1);
errs=zeros(Nh,3); itrs=zeros(Nh,1);
err=zeros(Nh,1);
tols=zeros(Nh,1); js=zeros(Nh,1);
for i=1:Nh
   js(i)=r(i);
   [u,v,p,xp,yp,xu,yu,xv,yv,x,y,~,itrs(i),fu,tols(i)] = StokesMAC(hs(i),f1,f2,utrue,vtrue,ptrue,40,'mg',r(i)-1); 
   ut=utrue(xu,yu);
   vt=vtrue(xv,yv);
   pt=ptrue(xp,yp);
   errs(i,1)=norm(abs(u(:)-ut(:)));
   errs(i,2)=norm(abs(v(:)-vt(:)));
   errs(i,3)=norm(abs(p(:)-pt(:)));
   err(i)=max(errs(i,1:2));
   Ns(i)=1/hs(i);
end
%% Plot Functions.
%
% xlabel('1/h'); ylabel('Iterations');
figure;
surf(xu,yu,u);
xlabel('x'); ylabel('y'); zlabel('u');
figure;
surf(xv,yv,v);
xlabel('x'); ylabel('y'); zlabel('v');
figure;
surf(xp,yp,p);
xlabel('x'); ylabel('y'); zlabel('p');

% figure;
% showrate3(Ns,errs(:,1),[],'-*','||u_I-u_h||_{0,2}', ...
%           Ns,errs(:,2),[],'-*','||v_I-v_h||_{0,2}', ...
%           Ns,errs(:,3),[],'-*','||p_I-p_h||_{0,2}');
%figure;
% showrate2(Ns,errs(:,1),[],'-*','||u_I-u_h||_{0,2}', ...
%           Ns,errs(:,2),[],'-*','||v_I-v_h||_{0,2}');
figure;
showrate(Ns,itrs,[],'-*','Iterations');

%% Save Results
%
% mgt=[hs,tols,itrs,err];
% dlmwrite('gausstests.txt', mgt , '&');
% mgt=[hs,js,tols,itrs,err];
% dlmwrite('tgtests.txt', mgt , '&');
mgt=[hs,js,tols,itrs,err];
dlmwrite('mgtests2.txt', mgt , '&');

% %% Plot Error
% %
% erru=abs(u-ut);
% errv=abs(v-vt);
% errp=abs(p-pt);
% figure;
% subplot(2,2,1);
% surf(xu,yu,erru);
% subplot(2,2,2);
% surf(xv,yv,errv);
% subplot(2,2,3);
% surf(xp,yp,errp);