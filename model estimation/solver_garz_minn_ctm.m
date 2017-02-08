function [erho,evel] = solver_garz_minn_ctm(tfinal,pos,D,Vel,cl,...
    cp,ca,h,w_min,w_max,rhom)

%==========================================================================
% This is a CTM like scheme for the second order models/not work for phast
% transition models. The model solved here is a generalized second order
% traffic model, which can be applied to other second order models by miner
% modifications.
% Feburary 11, 2014.
% Shimao Fan
% This is for test purpose, check weather the scheme works or not.
%*******************************************************************
%  U = [\rho,u], where rho is the density, and u is the velocity
%  Q = [rho;rho*w], where w is the property quantity
%  The inputs of this function involve
%  B = density data;
%  Vel = velocity data;
%  tfinal = length of simulation time
%  pos = gives the positions of sensor stations
%  cl, cp, ca = are the model parameters from fitting with data
%  h = is the size of spacial cell
%  rhom = the maximum density
%  w_min and w_max are bound of property quantity w, also from data fitting

%=========================================================================
%=========================================================================
% set up dt and dx
L = pos(end);     % length of road
N = ceil(L/h);    % number of grid
dx = L/N;         % grid size
%------------------------------------------
k = .95*h/w_max;     % approximate time step, w_max is maximun possible velocity
M = ceil(tfinal/k);         % number of time steps
dt = tfinal/M;              % here M is the number of time steps
lambda = dt/dx;      %size time cell over size of spacial cell
%------------------------------------------
nofdata = tfinal*120;       % this is number of data points/averaged in 30 seconds
inter = round(M/nofdata);   % how to pick the data
%------------------------------------------
% make a grid in time
t = 0:dt:tfinal;            % time vector, uniform
U=zeros(2,N);      % store the solution [density, velocity], i.e., [rho, u]
%------------------------------------------
% locate position of sensor 2/here we perform validation at sensor 2
position = ceil(pos(2:end-1)/dx);     % position of sensor 2
%------------------------------------------
% initial condition of density rho/kind of arbitrary
U(1,:)=10;     % initial density
U(2,:) = 100; % initial w, empty road velocity
%------------------------------------------
% The boundary conditoins based on inputs
%------------------------------------------
% find boundary condition by interpolation
% the cubic interpolation
%------------------------------------------
for i = 1:length(pos)   % number of sensor
    DensData(i,:) = intpol_bc(D(i,:),t);
    VelData(i,:) = intpol_bc(Vel(i,:),t);
end

% boundary condition
D_lbc = DensData(1,:); V_lbc = VelData(1,:);
D_rbc = DensData(end,:); V_rbc = VelData(end,:);





%///////////////////////
%  The main algorithm
%\\\\\\\\\\\\\\\\\\\\\\\
% intended to store the error data, measure the difference of model
% predition and data
e_rho = zeros(1,M);
e_vel = zeros(1,M);

Q = unknownQ(U);     % unknowns of in conserved form [rho, flux]
  
for n = 1:M      % loop over all time U = [density, w]; as unknowns
   % represent Q as function of U, the conservative unknowns
   %--------------------------------------------------------------
   %  Insert new boundary conditions on right hand side
   %  from velocity information to the w information, i.e., given u find w   
   rhol = D_lbc(n); rhor = D_rbc(n); % density of left and right bound
   vell = V_lbc(n); velr = V_rbc(n);
   
   %*******************************************************
   % solving inverse problem by applying newton's iteration
   %*******************************************************
   % given initial guess
   wl = 105;
   wr = wl; 
  %----------------------------------------------------  
   for i =1 : 2   %newton iteration/find w based on [rho u],density and velocity data
       wl = wl+(vell-Velocity(rhol,wl,cl,cp,ca,rhom))./...
           func_diff_v_w(rhol,wl,cl,cp,ca,rhom);
       wr = wr+(velr-Velocity(rhor,wr,cl,cp,ca,rhom))./...
          func_diff_v_w(rhor,wr,cl,cp,ca,rhom);       
   end
   %----------------------------------------------------   
   % perform projection, make sure w\in[w_min,w_max]
   wl = min(max(wl,w_min),w_max);
   wr = min(max(wr,w_min),w_max);
   %---------------------------------------------------
   % left and right boundary condition
   lbc = [D_lbc(n);wl];
   rbc = [D_rbc(n);wr];
   % right cell boundaries
   Up1=[U(:,2:end),rbc];    % right boundary
   % left cell boundaries
   Um1=[lbc,U(:,1:end-1)];  % left boundary 
   %---------------------------------------------------
   % uppdating rules based on sending and receiving
   Q=Q+lambda*(Num_Flux(Um1,U,cl,cp,ca,rhom)- ...
       Num_Flux(U,Up1,cl,cp,ca,rhom));
   %---------------------------------------------------
   %  find [rho,u] from conserved form, solution of next time step
   U(1,:)=Q(1,:);                     % traffic density
   U(2,:)=Q(2,:)./Q(1,:);             % traffic velocity
   %---------------------------------------------------

   rho = U(1,position);
   w = U(2,position);
   % find velocity from densityn rho and property w 
   vel = Velocity(rho,w,cl,cp,ca,rhom);
   % calculate the error
   data_rho = DensData(2:end-1,n)';
   data_vel = VelData(2:end-1,n)';
   
   e_rho(n) = norm(data_rho-rho,1);
   e_vel(n) = norm(data_vel-vel,1);

   %///////////////////////
   %   The plots part
   %///////////////////////
%     inter = 20;
%   if mod(n,inter) == 0
%    subplot(2,1,1)
%    plot(x,U(1,:),'r-','linewidth',3), hold on
%         plot(pos,DensData(n),'bs','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold on
%         plot(x(1),D_lbc(n),'bs','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold on
%         plot(x(end),D_rbc(n),'bs','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold off
%    title(sprintf('The traffic density at time t=%4.3f, w=%3.1f',n*dt,w),'fontsize',14)
%    axis([0,L,0,400])
%    ylabel('Density','fontsize',14)
% %    xlabel('Space','fontsize',14)
%    set(gca,'linewidth',2)
%    set(gca,'fontsize',14)
%    set(gca,'position',[.04 .55 .95 .4])
%    %----------------------------------------------------------------
%    subplot(2,1,2)
% %    plot(x,Velocity(U(1,:),U(2,:),x_lambda,x_p,rhomax),'b-','linewidth',2.5),hold on
%    plot(x,Velocity(U(1,:),U(2,:),x_lambda,x_p,rhomax),'b-','linewidth',3),hold on
% 
%        plot(pos,VelData(n),'bo','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold on
%        plot(x(1),V_lbc(n),'bo','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold on
%        plot(x(end),V_rbc(n),'bo','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold off
%    title(sprintf('The traffic velocity at t=%4.3f, w = %3.1f',n*dt,w),'fontsize',14)
%    axis([0,L,0,120])
%    ylabel('Velocity','fontsize',14)
%    xlabel('Space','fontsize',14)
%    set(gca,'linewidth',2)
%    set(gca,'fontsize',14)
%    res = 800;
%    set(gcf,'paperpositionmode','auto')
%    set(gcf,'position',[10  50 res res*.8])
%    set(gca,'position',[.04 .08 .95 .4])
%    drawnow
%    % Print figures
%    filename_save = sprintf('fig_garz_minn_day17_%06d',n/inter);
%    print(gcf,'-dpng',filename_save,'-r290')
%   end
end

% error calculation
start = 5*inter;


nn = length(pos)-2;    % number of sensors validated; the interior ones

erho = mean(e_rho(start:M))/(4*nn);
evel = mean(e_vel(start:M))/nn;

e_rho = e_rho(start:inter:end);
e_vel = e_vel(start:inter:end);

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
% partial derivative of w respect to Velocity function
%-------------------------------------------------------------
function y = func_diff_v_w(rho,w,cl,cp,ca,rhom)
        eps = 1.0e-04;
        dy = Velocity(rho,w+eps,cl,cp,ca,rhom) - ...
            Velocity(rho,w,cl,cp,ca,rhom);
        y = dy/eps;
%--------------------------------------------        
function  F = Num_Flux(Ul,Ur,cl,cp,ca,rhom) % numerical flux
   % find sending function
   s_flux = sending(Ul(1,:),Ul(2,:),cl,cp,ca,rhom);
   % receiving of flux
   r_flux = receiving(Ul,Ur,cl,cp,ca,rhom);
   % conservation of mass
   F(1,:) = min(s_flux,r_flux);   % take the minimum
   % conservation of momentum
   F(2,:) = Ul(2,:).*F(1,:);
%----------------------------------------------
% define the receiving function
function y = receiving(Ul,Ur,cl,cp,ca,rhom)
      w = Ul(2,:);
      p = polyval(cp,w);    % p parameter
      lambda = polyval(cl,w);  % \lambda parameter
%       alpha = polyval(ca,w);     % alpha parameter  
%       p = cp(1)*w.^2+cp(2)*w+cp(3);
%       lambda = cl(1)*w.^2+cl(2)*w+cl(3);
%       alpha = ca(1)*w.^2+ca(2)*w+ca(3);
      % determined by Ul
      rhoc = rho_critical(lambda,p,rhom);  % a vector of rhoc
      % velocity on the right side
      ur = Velocity(Ur(1,:),Ur(2,:),cl,cp,ca,rhom);
      % newton iteration is faster/less costy
      rho_middle = Riemann_G(Ul,Ur,cl,cp,ca,rhom); % middle state
      
      case1 = rho_middle<=rhoc;  % maximum 
      case2 = rho_middle>rhoc;   % maximum potential
      q_max = traffic_flux_smooth(rhoc,w,cl,cp,ca,rhom); % maximum flow rate
      y = case2.*rho_middle.*ur + case1.*q_max;  
%----------------------------------------------
% define sending functions
function y = sending(rho,w,cl,cp,ca,rhom)% U = [\rho, w] two state variables
%       p = cp(1)*w.^2+cp(2)*w+cp(3);
%       lambda = cl(1)*w.^2+cl(2)*w+cl(3);
      p = polyval(cp,w);    % p parameter
      lambda = polyval(cl,w);  % \lambda parameter
%       alpha = polyval(ca,w);     % alpha parameter
      rhoc = rho_critical(lambda,p,rhom);  % a vector of rhoc/3-params fd
      case1 = rho<=rhoc;  % maximum possible
      case2 = rho>rhoc;   % maximum
      q_max = traffic_flux_smooth(rhoc,w,cl,cp,ca,rhom); % maximum flow rate
      q = traffic_flux_smooth(rho,w,cl,cp,ca,rhom);
      y = case1.*q + case2.*q_max;
%-----------------------------------------------
% calculate the middle density
function rho = Riemann_G(Ul,Ur,cl,cp,ca,rhom) % calculate the riemann solution/intermediate state
 % Ul and Ur are initial data of Riemann problem
  u = Velocity(Ur(1,:),Ur(2,:),cl,cp,ca,rhom); % right velocity
  w = Ul(2,:);     % left propoty
  rho = Ul(1,:);   % not sure which one is a good pick
  for i =1 : 2   % something wrong with this root finding process
       rho = rho+(u-Velocity(rho,w,cl,cp,ca,rhom))./...
           DiffVel3(rho,w,cl,cp,ca,rhom);     
  end
%---------------------------------------------
function y = unknownQ(U)    % the conserved unknowns [rho, rho*w]
        y(1,:) = U(1,:);
        y(2,:) = U(1,:).*U(2,:);
%----------------------------------------------------         
function y = DiffVel3(x,w,cl,cp,ca,rhomax)% define derivative of velocity respect to density
%       p = cp(1)*w.^2+cp(2)*w+cp(3);
%       lambda = cl(1)*w.^2+cl(2)*w+cl(3);
%       alpha = ca(1)*w.^2+ca(2)*w+ca(3);
      p = polyval(cp,w);    % p parameter
      lambda = polyval(cl,w);  % \lambda parameter
      alpha = polyval(ca,w);     % alpha parameter
      a = sqrt(1+(p.*lambda).^2);
      c = ((x./rhomax)-p).*lambda;
      y = alpha.*(-lambda.*x.*c./(rhomax*sqrt(1+c.^2))-a+sqrt(1+c.^2))./(x.^2);
%---------------------------------   
function y = Velocity(x,w,x_lambda,x_p,x_alpha,rhomax) % velocity function 3 parameters
       % x is density;
   y = traffic_flux_smooth(x,w,x_lambda,x_p,x_alpha,rhomax)./x;   

function f = traffic_flux_smooth(x,w,cl,cp,ca,rhomax)% smooth fundamental diagram function
     % x_lambda and x_p are the coefficients of the fitting polynomials
      p = polyval(cp,w);    % p parameter
      lambda = polyval(cl,w);  % \lambda parameter
      alpha = polyval(ca,w);     % alpha parameter
%       p = cp(1)*w.^2+cp(2)*w+cp(3);
%       lambda = cl(1)*w.^2+cl(2)*w+cl(3);
%       alpha = ca(1)*w.^2+ca(2)*w+ca(3);
      
      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = ((x./rhomax)-p).*lambda;
      f = alpha.*(a+(b-a).*(x./rhomax)-sqrt(1+y.^2));      
%-------------------------------------------------------------   
function y = intpol_bc(data,t)   % of any order
    % data--- traffic data, traffic density
    % t --- time vector
    % assume picewise constant boundary conditions that came from real data
    tfinal = ceil(t(end)*60); % convert hour to minutes
    t = t*60;           % convert hour to minutes
    dt = 1/2;           % average avery 30 seconds, or half minute
    pl = 0:dt:tfinal;   % time interval we interested in for calculation
    n = length(pl);     % number of data we use
%     [param,res]= polyfit(pl,data,7);  % cl -- coefficient of lambda(w)
%     y = polyval(param,t);
    y = spline(pl,data(1:n),t);
%%%%%%%%%%%%%%%%%%%%%%End of
%%%%%%%%%%%%%%%%%%%%%%Code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

