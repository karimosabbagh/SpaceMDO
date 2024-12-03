function display_convergence(PB,x,v,w,obj)

q = get_q(PB,x);
q2 = q./PB.q_scale;
qmin = -v./(2*w.^2);
q2min = qmin./PB.q_scale;

[qmax,kmax] = max(abs(q)); % maximum inconsistency e_q

subplot(2,1,1); hold on;
if PB.no > 1
    h = findobj(gca,'Type','line'); 
    x = h.XData; 
    y = h.YData; 
    x = [x, PB.no];
    y = [y, qmax];
    cla(gca)
else
    x = [PB.no];
    y = [qmax];
end
    
plot(x,y,'-r');
ylabel('max inconsistency \epsilon_q');
xlabel('Outer loop iterations NO');
% set(gca,'yscale','log');
hold off;

subplot(2,1,2); hold on;
if PB.no > 1
    h = findobj(gca,'Type','line'); 
    x = h.XData; 
    y = h.YData; 
    x = [x, PB.no];
    y = [y, obj];
    cla(gca)
else
    x = [PB.no];
    y = [obj];
end
plot(x,y,'-b');
ylabel('Objective f(x)');
xlabel('Outer loop iterations NO');
hold off;

% subplot(2,1,1); hold on;
% x = [PB.no];
% y = [qmax];
% plot(x,y,'xr');
% ylabel('qmax');
% xlabel('no');
% hold off;
% 
% subplot(2,1,2); hold on;
% x = [PB.no];
% y = [obj];
% plot(x,y,'ob');
% xlabel('no');
% ylabel('obj');
% hold off;

drawnow

% % Capture the plot as an image 
% frame = getframe(gcf); 
% im = frame2im(frame); 
% [imind,cm] = rgb2ind(im,256); 
% % Write to the GIF File 
% if PB.no == 1 
%   imwrite(imind,cm,'progress.gif','gif', 'Loopcount',inf); 
% elseif PB.no > 1
%   imwrite(imind,cm,'progress.gif','gif','WriteMode','append'); 
% end 
