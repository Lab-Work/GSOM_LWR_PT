function [num_rho,num_vel,erho,evel] = solver_lwr_ngsim_godunov(tfinal,P,D,V,BD,x,vm,rhom,ref)

%==========================================================================
% Scalar model for Minneapolis data, the LWR model with data-fitted smooth
% 3-parameter fundamental diagram. It is a Godunov type solver.
%********************************************
% D,V --- is the density and velocity data
% lda, p, alpha --- fundamental diagram data, the equilibrium fd
% h --- grid size
% rhom, vm --- maximum density, and velocity
% Shimao Fan, Temple University
% modified on Jan 24 2013.
%==========================================================================
lda = P(1); p = P(2); alpha = P(3);
%-------------------------------------
% tfinal = 14;       % 14 minutes
h = x(2)-x(1);
m = length(x);     % number of grid points
N = length(x);
midd = round(N/2);
%--------------------------------
% time step
k = 1/36000; 
t0 = 0:k:(tfinal/60);
m1 = length(t0);
dt = 2*k/ref; 
% dt = k/ref;
t = 0:dt:(tfinal/60);
% boundary data/ only density for scalar model
lbc = spline(t0,BD(1:m1,1),t);
rbc = spline(t0,BD(1:m1,2),t);
% do projection
lbc(lbc>rhom) = rhom;
rbc(rbc>rhom) = rhom;
% interpol respect to time for all the data
for i = 1:size(D,2)   % loop respect space grids
    D_data(:,i) = spline(t0,D(1:m1,i),t);
    V_data(:,i) = spline(t0,V(1:m1,i),t);
end

D_data(D_data>rhom) = rhom;
% number of time step
M = ceil(tfinal/(60*dt));
% t = linspace(0,tfinal,N);
lambda = dt/h;

%///////////////////////
%  The Main Algorithm
%\\\\\\\\\\\\\\\\\\\\\\\

% start = 10*ceil(M/120);   % we not compare from t = 0,because initial condition 
erho = zeros(1,M);
evel = zeros(1,M);
U = D(1,:);   % initially uniform low density
U(U>=rhom) = rhom;
% calculate the critical density
rhoc = rho_critical(lda,p,rhom);

%///////////////////////
%  The Main Algorithm
%\\\\\\\\\\\\\\\\\\\\\\\


for j = 1:M   % loop with respect to time, each time step is 1/10 second
    
    ul = [lbc(j), U];   % length m+1 from 0 to m
    ur = [U,rbc(j)];    % length m+1 from 1 to m+1
    % specify the flux function on each piece of road
    f = @(rho) flux(rho,lda,p,alpha,rhom);

    fluxl = godunov(f,ul(1:m),ul(2:m+1),rhoc);
    fluxr = godunov(f,ur(1:m),ur(2:m+1),rhoc);
    % update to t+dt
    U = U + lambda * (fluxl - fluxr);
    % error calculation
    rho = U;
    f = flux(rho,lda,p,alpha,rhom);
    vel = f./rho;
    
    % store numerical result
    num_rho(j) = rho(midd);
    num_vel(j) = vel(midd);
    % error calculation
%     [e1,e2] = comput_error_nagim_time(D_data(j,:),...
%              V_data(j,:),rho,vel);
%     erho(j) = e1;
%     evel(j) = e2;

    rho_data = D_data(j,:);
    vel_data = V_data(j,:);
    rho_diff = rho-rho_data;
    vel_diff = vel-vel_data;
    erho(j) = norm(rho_diff(1:end),1);
    evel(j) = norm(vel_diff(1:end),1);
    
%     % do rescale
%     rho_diffs = (rho-rho_data)./rho_data;
%     vel_diffs = (vel-vel_data)./vel_data;
%     erhos(j) = norm(rho_diffs,1);
%     evels(j) = norm(vel_diffs,1);
    %==============================
    %       plotting parts
%     %==============================  
%    if mod(j,100) == 0
%    subplot(2,1,1),
%      plot([px(1:2:end)';px(end)],rho,'m-','linewidth',3), hold on
%      plot(pos,DensData(j),'bs','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold on
%         plot(px(1),lbc(j),'bs','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold on
%         plot(px(end),rbc(j),'bs','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold off
% 
% %      plot(x,DD(j,:),'b-','linewidth',2.5), hold off
%      axis([0 len 0 rhom])
%      legend('Num','Emp')
%      set(gca,'linewidth',2)
%      set(gca,'fontsize',14)
%      title(sprintf('traffic density t =%3.0fmin  ',j*dt*60),'fontsize',14)
%      legend('numerical sol','empirical data')
%      set(gca,'xtick',[])
%      ylabel('Density (# of vehicles/km)','fontsize',14)
%      xlabel('km','fontsize',14)
%      set(gca,'fontsize',14)
% %          ===============
% %           Make Movie
% %          ===============
% %          plot of errors
%   subplot(2,1,2),
%      plot([px(1:2:end)';px(end)],vel,'m-','linewidth',3), hold on
%         plot(pos,VelData(j),'bo','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold on
%        plot(px(1),V_lbc(j),'bo','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold on
%        plot(px(end),V_rbc(j),'bo','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold off
% 
% %      plot(t(1:50:end),e(1:50:end),'k--.','linewidth',1)
%      set(gca,'linewidth',2)
%      set(gca,'fontsize',14)
%    title(sprintf('The traffic velocity at time t=%3.0f',j*dt*60))
%    axis([0,len,0,120])
%    ylabel('Velocity','fontsize',14)
%    xlabel('Space','fontsize',14)     
%    set(gca,'fontsize',14)
%      res = [971 600];
%      set(gcf,'paperpositionmode','auto')
%      set(gcf,'position',[10  50 res(1) res(2)])
%      drawnow
% %      filename_save = sprintf('lwr_day38_%06d',j/10);
% %      print(gcf,'-dpng',filename_save,'-r200')        
%     end   % the end of 'if'
%--------------------------------------------------------------------------
end

erho = erho/ref;
evel = evel/ref;

% %\\\
% erhos = erhos/ref;
% evels = evels/ref;
%fprintf('\n')
%==========================================================================     
%==========================================================================     
%==========================================================================     
%==========================================================================     
%==========================================================================     
%==========================================================================
    
% function y = intpol_bc(data,t)
%     % data--- traffic data, traffic density
%     % t --- time vector
%     % assume picewise constant boundary conditions that came from real data
%     tfinal = t(end)*60; % convert hour to minutes
%     t = t*60;           % convert hour to minutes
%     dt = 1/2;           % average avery 30 seconds, or half minute
%     pl = 0:dt:tfinal;   % time interval we interested in for calculation
%     n = length(pl);     % number of data we use
%     y = spline(pl,data(1:n),t);
    
%------------------------------------------------------------
function z = godunov(fct,a,b,c)  % godunov flux better for nonlinear flux
% fct is the flux function
% where c is the critical points
case1 = a <= b;
case2 = a > c & b < c;
case3 = not(case1 | case2);
z = case1 .* min(feval(fct,a),feval(fct,b)) + ...
    case2 .*feval(fct,c) + case3 .* max(feval(fct,a),feval(fct,b));
%---------------------------------------------------------
% determine the critical density, for any given smooth curve
function y = rho_critical(lambda,p,rhom)
%       lambda = polyval(p1,x);
%       p = polyval(p2,x);
      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = rhom*((b-a)./(lambda.*(sqrt(lambda.^2-(b-a).^2)))+p);
%---------------------------------------------------------
function y = flux(rho,lambda,p,alpha,rhom) % velocity function/equilibrium      
%       lambda = polyval(p1,x);
%       p = polyval(p2,x);
%       alpha = polyval(p3,x);
%       whos
      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = ((rho./rhom)-p).*lambda;
      y = alpha.*(a+(b-a).*(rho./rhom)-sqrt(1+y.^2));
%=============================end of the code==============================
