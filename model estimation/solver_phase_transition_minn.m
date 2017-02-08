function [erho,evel] = solver_phase_transition_minn(tfinal,pos,D,Vel,u_r,...
   rho_r,h,w_max,q_min,q_max,rhom)
%==========================================================================
% solver for phase transition model, modified version.
%==========================================================================
% close all;
% some parameters
% pick the middle curve/ NGSIM data
% lda = P(1);
% p = P(2);
% alpha = P(3);
% rhovec = .1:.1:rhom;           % traffic density vector       
%------------------------------------------------
% determine dx and dt
L = pos(end);
% tfinal = 2.0;     % final time/hour
N = ceil(L/h);    % number of grid
dx = L/N;         % final grid size
%------------------------------------------
% k = .9*h/sl;   % approximate time step
k = .95*h/w_max;     % approximate time step, vm is maximun velocity
M = ceil(tfinal/k);         % number of time steps
dt = tfinal/M;              % here M is the number of time steps
lambda = dt/dx;
%------------------------------------------
nofdata = tfinal*120;   % this is number of data points
inter = round(M/nofdata);   % how to pick the data
%------------------------------------------
%------------------------------------------
% make a grid
%--------------
x = transpose(dx/2:dx:L-dx/2)'; % stagged grids
t = 0:dt:tfinal;                % time vector, uniform
U=zeros(2,N);      % the solution [density, velocity], i.e., [rho, u]
%------------------------------------------
% locate position of sensor 2
position = ceil(pos(2:end-1)/dx);     % position of sensor 2
%------------------------------------------
% initial condition of density rho
%-----------------------------------
U(1,:)=10;        % initial density
U(2,:) = 5*1.0e+03; % initial capacity
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

%*************************************************************************
%*******************  ----------------------  ****************************
%                       The main algorithm
%*******************  ----------------------  ****************************
%*************************************************************************

e_rho = zeros(1,M);
e_vel = zeros(1,M);
% w_min = 55;
% w_max = 80;
Q = unknownQ(U);     % unknowns of in conserved form [rho, flux]
% q_min = 2000;

for n = 1:M      % loop over all time U = [density, w]; as unknowns
   %  Insert new boundary conditions on right hand side
   rhol = D_lbc(n); rhor = D_rbc(n); % density of left and right bound
   ul = V_lbc(n); ur = V_rbc(n);
   
   % apply bisection root finder
   ql = func_inverse_v_pase_trans(rhol,ul,q_min,q_max,u_r,rho_r,rhom);
   qr = func_inverse_v_pase_trans(rhor,ur,q_min,q_max,u_r,rho_r,rhom);
   
%    ql = min(max(ql,q_min),q_max);
%    qr = min(max(qr,q_min),q_max);

   lbc = [D_lbc(n);ql];
   rbc = [D_rbc(n);qr];
   Up1=[U(:,2:end),rbc];             % right boundary
   Um1=[lbc,U(:,1:end-1)];           % left boundary 

   %---------------------------------------------------
   Q=Q+lambda*(HLL(Um1,U,u_r,rho_r,rhom)- ...
       HLL(U,Up1,u_r,rho_r,rhom));   
   %----------------------------------------------------------------------
   %  find [rho,u] from conserved form, solution of next time step
   %-----------------------------------------------------------------
   U(1,:)=Q(1,:);                     % traffic density
   U(2,:)=Q(2,:)./Q(1,:);             % traffic velocity
   %----------------------------------------
   rho = U(1,position);
   q = U(2,position);
%    beta = min(max(beta,beta_min),beta_max);

%    w = projection(w,w_min,w_max);
   vel = func_u_phase_transition(rho,q,u_r,rho_r,rhom);
   %----------------------------------------

   % error calculation
   data_rho = DensData(2:end-1,n)';
   data_vel = VelData(2:end-1,n)';
   diff_rho = data_rho-rho;
   diff_vel = data_vel-vel;
   
   e_rho(n) = norm(diff_rho,1);
   e_vel(n) = norm(diff_vel,1);
   %////////////////////
   %   The plots part
   %\\\\\\\\\\\\\\\\\\\\
%   if mod(n,10) == 0
%    subplot(2,1,1)
%    plot(x,U(1,:),'m-','linewidth',2.5), hold on
%    plot(x,D(n,:),'b-','linewidth',2.5), hold off
%    title(sprintf('The traffic density at time t=%4.3f',n*dt))
%    axis([0,x(end),0,600])
%    ylabel('Density','fontsize',14)
%    xlabel('Space','fontsize',14)
%    set(gca,'linewidth',2)
%    set(gca,'fontsize',14)
% %    %----------------------------------------------------------------
%    subplot(2,1,2)
%    plot(x,vel,'m-','linewidth',2.5),hold on
%    plot(x,V(n,:),'b-','linewidth',2.5), hold off
% %    title(sprintf('The traffic velocity at time t=%4.3f',n*dt))
%    axis([0,L,0,100])
%    ylabel('Velocity','fontsize',14)
%    xlabel('Space','fontsize',14)
%    set(gca,'linewidth',2)
%    set(gca,'fontsize',14)
%    res = [971 600];
%    set(gcf,'paperpositionmode','auto')
%    set(gcf,'position',[10 50 res(1) res(2)])
%    drawnow
% %    Print figures
% %    filename_save = sprintf('AwRscal_%06d',n/40);
% %    print(gcf,'-dpng',filename_save,'-r200')
%   end

%------------------------------------------------
% Calculate error 
%    rho_diff = rho-D_data(n,:);
%    vel_diff = vel-V_data(n,:);
% %     
%    e1 = norm(rho_diff,1);
%    e2 = norm(vel_diff,1);
%          
%    erho(n) = e1;
%    evel(n) = e2;
%    n
end

%---------------------------------
% error calculation
% start = 10*ceil(M/120);
start = 5*inter;
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
function y = func_diff_v_beta(rho,beta,ct,sl,rhom)
        eps = 1.0e-09;
        dy = Velocity(rho,beta+eps,ct,sl,rhom) - ...
            Velocity(rho,beta,ct,sl,rhom);
        y = dy/eps; 
%---------------------------------------------------------------
        
function F = HLL(Ul,Ur,u_r,rho_r,rhom) % the HLL Riemann solver for fds of 3 parameters
        % eigenvalues:
        eigns1 = lambda(Ul,u_r,rho_r,rhom);  % left bound
        eigns2 = lambda(Ur,u_r,rho_r,rhom);  % right bound
        % another way to define c1 and c2
        c1 = minimum(eigns1(1,:),eigns2(1,:));
        c2 = maximum(eigns1(2,:),eigns2(2,:));
        case1 = c1>=0;
        case2 = c1<0 & c2>=0;
        case3 = c2<0;
        Ql = unknownQ(Ul);
        Qr = unknownQ(Ur);
        % left flux
        Fl = flux(Ul,u_r,rho_r,rhom);
        % right flux
        Fr = flux(Ur,u_r,rho_r,rhom);
        % intermediate flux
        F_hll(1,:) = (c2.*Fl(1,:)-c1.*Fr(1,:)+(Qr(1,:)-Ql(1,:)).*c1.*c2)./...
            (c2-c1);
        F_hll(2,:) = (c2.*Fl(2,:)-c1.*Fr(2,:)+(Qr(2,:)-Ql(2,:)).*c1.*c2)./...
            (c2-c1);
        F(1,:) = case1.*Fl(1,:)+case2.*F_hll(1,:)+case3.*Fr(1,:);
        F(2,:) = case1.*Fl(2,:)+case2.*F_hll(2,:)+case3.*Fr(2,:);
%--------------------------------------------------------------------------   
function y = lambda(U,u_r,rho_r,rhom)  % U = (density, w), calculate eigenvalues
       v = func_u_phase_transition(U(1,:),U(2,:),u_r,rho_r,rhom); % velocity
%        diffv=diff_vel_pt(U(1,:),U(2,:),u_r,rhom);
       diffv = DiffVel(U(1,:),U(2,:),u_r,rho_r,rhom);
       y(1,:) = v + U(1,:).*diffv;   % modified here
       y(2,:) = v;
%----------------------------------------------
function y = flux(U,u_r,rho_r,rhom)  % U = [rho,beta], the numerical flux
        v = func_u_phase_transition(U(1,:),U(2,:),u_r,rho_r,rhom);    % velocity
        y(1,:) = U(1,:).*v;
        y(2,:) = U(1,:).*U(2,:).*v;
%---------------------------------------------
function y = unknownQ(U)    % the conserved unknowns [rho, rho*beta]
        y(1,:) = U(1,:);
        y(2,:) = U(1,:).*U(2,:);
%----------------------------------------------------   
function y = DiffVel(rho,q,u_r,rho_r,rhom)% define derivative of velocity
      eps = 1.0e-05;
      y =(func_u_phase_transition(rho+eps,q,u_r,rho_r,rhom)-...
          func_u_phase_transition(rho,q,u_r,rho_r,rhom))/eps;
% -------------------------------------------------------      
% function y = Velocity(rho,beta,ct,sl,rhom) % velocity function 3 parameters
%        rhoc = polyval(ct,beta);
%        y = traffic_flux_trig(rho,sl,rhoc,rhom)./rho; 
%---------------------------------------------------------------- 
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

