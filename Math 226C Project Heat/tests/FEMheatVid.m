%% Finite Element Tests
%
%   Written by Ted Kwan for Math 226C Project Heat
%
%   This script runs the tests for the finite element method,
%   and saves the video to plot.

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

dt=0.01; t0=0; tf=2; h=0.02; % Space Error
[node,elem] = squaremesh([0,1,0,1],h);
[ut,node,elem,t,M,A,At1]=FiniteElemHeatVideo(node,elem,f,u0,t0,tf,dt,utrue);
% mv = VideoWriter('HeatUniformRes.avi');
% mv.FrameRate=4; open(mv);
fig=figure;
for k = 1:5:length(ut)
    %%% Plot Figures
    %
    showresult(node,elem,ut{k},[-9,20]); axis([0,1,0,1,-0.01,0.1]);
    legend(['t=' num2str(t(k))],'LOCATION','best');
    pause(0.01);
    %%% Uncomment to save movie
% 	frame = getframe(fig); % Capture t spot.
%     writeVideo(mv,frame); % Save.
end
% %%% Write Movie.
% %
% close(mv);