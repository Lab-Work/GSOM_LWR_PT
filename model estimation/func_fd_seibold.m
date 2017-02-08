function f = func_fd_seibold(rho,a,b,v0,rho0,f0,rhom)

% here, a -- controal the curvature
%       b -- maximum point, or critical density
% [u0,v0,rho0] are initial flux, velocity, and density, which are fd
% parameters.
% Idea of Dr. Benjamin Seibold,

% May 26, 2013
% Shimao Fan, Temple University
%---------------------------------------------------------


if nargin<1
% initial flux
f0 = 5000;   
% free flow density
rho0 = 56; 
% initial velocity
v0 = 70;
rhom = 800;  % the maximum density
rho = rho0:rhom;
a = 10;
b = 1.5*rho0;
end

drho = rhom-rho0;

w = @(rho) -atan((rho-b)/a);
w0 = w(rho0);
%-------------------------------------------
% transform w to v
I = int_actan(rhom,a,b)-int_actan(rho0,a,b);
k = (v0*drho+f0)/(w0*drho-I); % slop
inters = v0-k*w0;   % y intersection
%-------------------------------------------
 % derivative of curve
    % scaling and shifting so that v(0)=1 and \int v(x) dx = -u0
% v = k*w + inters;
%-------------------------------------------
f1 = @(rho) k*int_actan(rho,a,b)+inters*rho;
C = f0-f1(rho0);   %value of constant term
f = f1(rho)+C;  

%=====================================
% define a curve family
%-------------------------------------
%\\\\\\\\\\\\\
%plotting part
%/////////////
% figure(1)
% clf
% close all;
% % vertices of the triangle
% vert_x = [rho0 rhom rhom];
% vert_y = [u0 func_ub(rhom,v0,rho0,u0) 0];
% %----------------------------
% % construct internal curves
% a = 10.^(linspace(-4,3,10));
% b = linspace(rho0,rhom,10);
% 
% [A,B] = meshgrid(a,b);
% for i = 1:length(a)
%     for j = 1:length(b)
%         plot()
%     end
% end
% f_ub = func_ub(rho,v0,rho0,u0);   % up bound of function
% f_lb = func_lb(rho,rho0,u0,rhom); % low bound of function
% 
% plot(rho,f_ub,'-.','color',[0 0 0.7],'linewidth',4), hold on
% plot(rho,f_lb,'-.','color',[0 0 0.7],'linewidth',4), hold on
% plot(rhom+0*(0:vert_y(2)),0:vert_y(2),'-.','color',[0 0 0.7],'linewidth',4), hold on
% plot(rho,f,'k-','linewidth',2), hold on
% plot(vert_x,vert_y,'.','color',[0.8 0 0],'markersize',30), hold on
% % plot(rho0:rhom,vert_y(2)+0*(rho0:rhom),'k--','linewidth',4), hold on
% axis([rho0 1.1*rhom 0 1.1*max(f_ub)])
% xlabel('density \rho','fontsize',14)
% ylabel('flow rate','fontsize',14)
% set(gca,'linewidth',2,'fontsize',14)
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% 
% % text(rho0+10,-1000,'\rho_{f}','fontsize',14,'horizontalalign','center')
% text(rhom,-1500,'(\rho_{max},0)','fontsize',14,'horizontalalign','center')
% text(rho0-50,u0,'(\rho_{f},q_{f})','fontsize',14)
% text(rhom,1.04*vert_y(2),'(\rho_{max},v_{f}(\rho_{max}-\rho_{f})+q_{f})','fontsize',14,'horizontalalign','center')
% 
% % text(rho0-15,0,'0','fontsize',14)
% 
% res = 600;
% set(gcf,'paperpositionmode','auto')
% set(gcf,'position',[10 50 res*1.15 res])
% set(gca,'position',[.06 .08 .92 .88])
% title('flux function')


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%///////////////////////////////////////////////
function y = func_ub(rho,v0,rho0,f0)
y = v0*(rho-rho0)+f0;

function y = func_lb(rho,rho0,f0,rhom)

s = -f0/(rhom-rho0);  % negative slop
y = s*(rho-rhom);

function y = int_actan(x,a,b)
% y = -a*((x2-b)/a*atan((x2-b)/a)-(x1-b)/a*atan((x1-b)/a)+0.5*...
%     log((a^2+(x1-b)^2)/(a^2+(x2-b)^2)));
y = -a*((x-b)/a.*atan((x-b)/a)-0.5*log(1+(x-b).^2/(a^2)));
