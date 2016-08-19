figure;
plot(Ns,timen,Ns,timetg,Ns,timefas);
xlabel('\frac{1}{h}'); ylabel('times');
title('CPU Time Used for Different Methods');
axis tight; grid on;
h_legend=legend('Newton','Two Grid','FAS','LOCATION','Best');
set(h_legend,'FontSize',12);