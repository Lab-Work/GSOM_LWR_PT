% function [rho_num,vel_num,erho,evel] = solver_lwrq_minn_godunov(pos,D,Vel,vm,wm,h,rhom)
function [erho,evel] = solver_lwrq_minn_godunov(tfinal,pos,D,Vel,vm,wm,h,rhom)

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
% generate dx, and dt
% len = 1.224;          % length of the study area
len = pos(end);
% tfinal = 1.0;           % time inverval 1 hour
% m = ceil(len./h);     % number of grid points
approxk = .95*h/wm;    % CFL condition
N = ceil(tfinal/approxk);
k = tfinal/N; dt = k;   % time step
lambda = k/h;
t = 0:dt:tfinal;        % t vector for boundary condition
%------------------------------------------
nofdata = tfinal*120;   % this is number of data points
inter = round(N/nofdata);   % how to pick the data
%------------------------------------------
%------------------------------------------
x = transpose(h/2:h:len-h/2)';
m = length(x);
position = ceil(pos(2:end-1)/h);     % position of sensor 2
%--------------------------------------------
% calculate critical density
rhoc = rhom/2;
%--------------------------------------------
% interpolate with respect to time
for i = 1:length(pos)   % number of sensor
    DensData(i,:) = intpol_bc(D(i,:),t);
    VelData(i,:) = intpol_bc(Vel(i,:),t);
end

% the boundary condition;
lbc = DensData(1,:);   %V_lbc = VelData(1,:);
rbc = DensData(end,:); %V_rbc = VelData(end,:);


lbc(lbc>rhom) = rhom;
rbc(rbc>rhom) = rhom;

DensData(DensData>rhom) = rhom;
%///////////////////////
%  The Main Algorithm
%\\\\\\\\\\\\\\\\\\\\\\\

% start = 10*ceil(N/240);   % we not compare from t = 0,because initial condition 
start = 5*inter;   % we not compare from t = 0,because initial condition 
e_rho = zeros(1,N);
e_vel = zeros(1,N);
% rho_num = zeros(1,N);
% vel_num = zeros(1,N);

U = zeros(1,m);   % initial condition

for j = 1:N    % loop with respect to time, each time step is 1/10 second
    
    ul = [lbc(j), U];   % length m+1 from 0 to m
    ur = [U,rbc(j)];    % length m+1 from 1 to m+1
    % specify the flux function on each piece of road
    f = @(rho) flux(rho,vm,rhom);

    fluxl = godunov(f,ul(1:m),ul(2:m+1),rhoc);
    fluxr = godunov(f,ur(1:m),ur(2:m+1),rhoc);
    % update to t+dt
    U = U + lambda * (fluxl - fluxr);
    % error calculation
    rho = U(position);
    f = flux(rho,vm,rhom);
    vel = f./rho;
    
%     rho_num(j) = rho;
%     vel_num(j) = vel;
       
    data_rho = DensData(2:end-1,j)';
    data_vel = VelData(2:end-1,j)';
    diff_rho = data_rho-rho;
    diff_vel = data_vel-vel;
   
    e_rho(j) = norm(diff_rho,1);
    e_vel(j) = norm(diff_vel,1);
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
 
% erho = mean(abs(DensData(start:N)-NumDensity(start:N)))/rhom;
% evel = mean(abs(VelData(start:N)-NumVel(start:N)))/vm;
% sub functions used for this code

nn = length(pos)-2;    % number of sensors validated; the interior ones

erho = mean(e_rho(start:N))/nn;
evel = mean(e_vel(start:N))/nn;


e_rho = e_rho(start:inter:end);
e_vel = e_vel(start:inter:end);
%fprintf('\n')
%==========================================================================     
%==========================================================================     
%==========================================================================     
%==========================================================================     
%==========================================================================     
%==========================================================================
    
function y = intpol_bc(data,t)
    % data--- traffic data, traffic density
    % t --- time vector
    % assume picewise constant boundary conditions that came from real data
    tfinal = t(end)*60; % convert hour to minutes
    t = t*60;           % convert hour to minutes
    dt = 1/2;           % average avery 30 seconds, or half minute
    pl = 0:dt:tfinal;   % time interval we interested in for calculation
    n = length(pl);     % number of data we use
    y = spline(pl,data(1:n),t);
    
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

      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = rhom*((b-a)./(lambda.*(sqrt(lambda.^2-(b-a).^2)))+p);
%---------------------------------------------------------
function y = flux(rho,vm,rhom) % velocity function/equilibrium      
    y = vm*rho.*(1-rho/rhom);
%=============================end of the code==============================
