function [erho,evel] = solver_garz_minn(tfinal,pos,D,Vel,cl,...
    cp,ca,h,w_min,w_max,rhom)

%==========================================================================
%  Simulate AW-Rascal model by Lax-friendrish and HLL approximate Riemann 
%  solver method, In this code, we take the data into classical AR model,
%  that is, the initial and boundary condition are from real traffic data.
%*******************************************************************
%  U = [\rho,u], where rho is the density, and u is the velocity
%  Q = [rho;rho*(u+P(rho))], where P(rho) is pressure function
%  B = density data;
%  Vel = velocity data;
%  [a,b] is the time intervel, used to select proper traffic data
%  x_lambda, x_p and x_alpha are parameter of fitting results
%  Mw is the matrix store corresponding data of "rho-w-v"
%  h = size of grid point
%  rhomax = the maximum density
%*******************************************************************
%  Shimao Fan, Jan 24, 2012.
%  @ Temple University, Math Department.
%  Modified on Jan 24 2013
%=========================================================================
% rhovec = .1:.1:rhomax;           % traffic density vector
%--------------------------------------------
% set up dt and dx
L = pos(end);
% tfinal = 1.0;     % final time/hour
% tfinal = 2.0;
N = ceil(L/h);    % number of grid
dx = L/N;         % final grid size
%------------------------------------------
% k = .5*h/max(Vel(1,a:b));   % approximate time step
k = .95*h/w_max;     % approximate time step, vm is maximun velocity
M = ceil(tfinal/k);         % number of time steps
dt = tfinal/M;              % here M is the number of time steps
lambda = dt/dx;
%------------------------------------------
nofdata = tfinal*120;   % this is number of data points
inter = round(M/nofdata);   % how to pick the data
%------------------------------------------
% make a grid
%--------------
x = transpose(dx/2:dx:L-dx/2)'; % stagged grids
t = 0:dt:tfinal;                % time vector, uniform
U=zeros(2,N);      % the solution [density, velocity], i.e., [rho, u]
%------------------------------------------
% locate position of sensor 2
% pos = .616;
% position = ceil(pos/dx);     % position of sensor 2
position = ceil(pos(2:end-1)/dx);     % position of sensor 2
%------------------------------------------
% initial condition of density rho
%-----------------------------------
U(1,:)=10;     % initial density
U(2,:) = 100; % initial w, empty road velocity
%------------------------------------------
% The boundary conditoins
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
% the dt will change each time step corresponding to the CFL condition,
% this will increase the difficulty to install boundary conditions
e_rho = zeros(1,M);
e_vel = zeros(1,M);
% rho_num = zeros(1,M);
% vel_num = zeros(1,M);

Q = unknownQ(U);     % unknowns of in conserved form [rho, flux]
  
for n = 1:M      % loop over all time U = [density, w]; as unknowns
   % represent Q as function of U, the conservative unknowns
   %--------------------------------------------------------------
   %  Insert new boundary conditions on right hand side
   %  from velocity information to the w information, i.e., given u find w   
   rhol = D_lbc(n); rhor = D_rbc(n); % density of left and right bound
   vell = V_lbc(n); velr = V_rbc(n);
   % given initial guess
   % apply newton's iteration
%    wl = 1.1*w_max;
   wl = 105;
%    wl = 94;
   wr = wl; 
%  %----------------------------------------------------  
   for i =1 : 2   % something wrong with this root finding process
       wl = wl+(vell-Velocity3(rhol,wl,cl,cp,ca,rhom))./...
           func_diff_v_w(rhol,wl,cl,cp,ca,rhom);
       wr = wr+(velr-Velocity3(rhor,wr,cl,cp,ca,rhom))./...
          func_diff_v_w(rhor,wr,cl,cp,ca,rhom);       
   end
   %----------------------------------------------------   
      % apply bisection root finder
%    wl = func_inverse_v(rhol,vell,w_min,w_max,cl,cp,ca,rhom);
%    wr = func_inverse_v(rhor,velr,w_min,w_max,cl,cp,ca,rhom);
   
   % perform projection
   wl = min(max(wl,w_min),w_max);
   wr = min(max(wr,w_min),w_max);
   %---------------------------------------------------
   lbc = [D_lbc(n);wl];
   rbc = [D_rbc(n);wr];
   Up1=[U(:,2:end),rbc];    % right boundary
   Um1=[lbc,U(:,1:end-1)];  % left boundary 
   %---------------------------------------------------
   % HLL Riemann solver
   Q=Q+lambda*(HLL3(Um1,U,cl,cp,ca,rhom)- ...
       HLL3(U,Up1,cl,cp,ca,rhom));
   %---------------------------------------------------
   %  find [rho,u] from conserved form, solution of next time step
   U(1,:)=Q(1,:);                     % traffic density
   U(2,:)=Q(2,:)./Q(1,:);             % traffic velocity
   %---------------------------------------------------
%    U(:,1) = lbc;
%    U(:,end) = rbc;
   rho = U(1,position);
   w = U(2,position);
   % perform projection process  
   vel = Velocity3(rho,w,cl,cp,ca,rhom);
   
   data_rho = DensData(2:end-1,n)';
   data_vel = VelData(2:end-1,n)';
   
   e_rho(n) = norm(data_rho-rho,1);
   e_vel(n) = norm(data_vel-vel,1);
   
   
%    rho_num(n) = rho;
%    vel_num(n) = vel;
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
% start = 10*ceil(M/240);


% erho = mean(abs(DensData(start:M)-NumDens(start:M)))/400;
% evel = mean(abs(VelData(start:M)-NumVel(start:M)))/120;
% erho = mean(abs(DensData(start:M)-NumDens(start:M)))/rhom;
% evel = mean(abs(VelData(start:M)-NumVel(start:M)))/vm;

nn = length(pos)-2;    % number of sensors validated; the interior ones

erho = mean(e_rho(start:M))/nn;
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
        dy = Velocity3(rho,w+eps,cl,cp,ca,rhom) - ...
            Velocity3(rho,w,cl,cp,ca,rhom);
        y = dy/eps;
 %------------------------------------------------------------------       
function F = HLL3(Ul,Ur,x_lambda,x_p,x_alpha,rhomax) % the HLL Riemann solver for fds of 3 parameters
        % eigenvalues:
        eigns1 = lambda3(Ul,x_lambda,x_p,x_alpha,rhomax);  % left bound
        eigns2 = lambda3(Ur,x_lambda,x_p,x_alpha,rhomax);  % right bound
        % another way to define c1 and c2
        c1 = minimum(eigns1(1,:),eigns2(1,:));
        c2 = maximum(eigns1(2,:),eigns2(2,:));
        case1 = c1>=0;
        case2 = c1<0 & c2>=0;
        case3 = c2<0;
        Ql = unknownQ(Ul);
        Qr = unknownQ(Ur);
        % left flux
        Fl = flux3(Ul,x_lambda,x_p,x_alpha,rhomax);
        % right flux
        Fr = flux3(Ur,x_lambda,x_p,x_alpha,rhomax);
        % intermediate flux
        F_hll(1,:) = (c2.*Fl(1,:)-c1.*Fr(1,:)+(Qr(1,:)-Ql(1,:)).*c1.*c2)./...
            (c2-c1);
        F_hll(2,:) = (c2.*Fl(2,:)-c1.*Fr(2,:)+(Qr(2,:)-Ql(2,:)).*c1.*c2)./...
            (c2-c1);
        F(1,:) = case1.*Fl(1,:)+case2.*F_hll(1,:)+case3.*Fr(1,:);
        F(2,:) = case1.*Fl(2,:)+case2.*F_hll(2,:)+case3.*Fr(2,:);
%--------------------------------------------------------------------------   
function y = lambda3(U,x_lambda,x_p,x_alpha,rhomax)  % U = (density, w), calculate eigenvalues
       v = Velocity3(U(1,:),U(2,:),x_lambda,x_p,x_alpha,rhomax); % velocity
       diffv = DiffVel3(U(1,:),U(2,:),x_lambda,x_p,x_alpha,rhomax);
       y(1,:) = v + U(1,:).*diffv;   % modified here
       y(2,:) = v;
%----------------------------------------------        
function y = flux3(U,x_lambda,x_p,x_alpha,rhomax)  % U = [rho,w], the numerical flux
        v = Velocity3(U(1,:),U(2,:),x_lambda,x_p,x_alpha,rhomax);    % velocity
        y(1,:) = U(1,:).*v;
        y(2,:) = U(1,:).*U(2,:).*v;
%---------------------------------------------
function y = unknownQ(U)    % the conserved unknowns [rho, rho*w]
        y(1,:) = U(1,:);
        y(2,:) = U(1,:).*U(2,:);
%----------------------------------------------------         
function y = DiffVel3(x,w,cl,cp,ca,rhomax)% define derivative of velocity
%       p = polyval(cp,w); 
%       lambda = polyval(cl,w);  % input data
%       alpha = polyval(ca,w);     % alpha parameter
      p = cp(1)*w.^2+cp(2)*w+cp(3);
      lambda = cl(1)*w.^2+cl(2)*w+cl(3);
      alpha = ca(1)*w.^2+ca(2)*w+ca(3);
      a = sqrt(1+(p.*lambda).^2);
      c = ((x./rhomax)-p).*lambda;
      y = alpha.*(-lambda.*x.*c./(rhomax*sqrt(1+c.^2))-a+sqrt(1+c.^2))./(x.^2);
%---------------------------------   
function y = Velocity3(x,w,x_lambda,x_p,x_alpha,rhomax) % velocity function 3 parameters
       % x is density;
   y = traffic_flux_smooth3(x,w,x_lambda,x_p,x_alpha,rhomax)./x;   

function f = traffic_flux_smooth3(x,w,cl,cp,ca,rhomax)% smooth flux function
     % x_lambda and x_p are the coefficients of the fitting polynomials
%       p = polyval(cp,w);    % p parameter
%       lambda = polyval(cl,w);  % \lambda parameter
%       alpha = polyval(ca,w);     % alpha parameter
      p = cp(1)*w.^2+cp(2)*w+cp(3);
      lambda = cl(1)*w.^2+cl(2)*w+cl(3);
      alpha = ca(1)*w.^2+ca(2)*w+ca(3);
      
      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = ((x./rhomax)-p).*lambda;
      f = alpha.*(a+(b-a).*(x./rhomax)-sqrt(1+y.^2));      
%-------------------------------------------------------------
function y = minimum(a,b) % here a and b could be vectors
     temp1 = a<=b;
     temp2 = a>b;
     y = temp1.*a+temp2.*b;    % where y is the minimum
  
function y = maximum(a,b) % so as above function
     temp1 = a<=b;
     temp2 = a>b;
     y = temp1.*b+temp2.*a;
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
%%%%%%%%%%%%%%%%%%%%%%End of Code%%%%%%%%%%

