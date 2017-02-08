function [num_rho,num_vel,erho,evel] = solver_garz_ngsim_collapsed_ctm(tfinal,D,V,BD,BV,ca,cb,...
   rho0,v0,f0,uf,rhof,x,ref,q_min,q_midd,q_max,rhom)
proj = 1;
%==========================================================================
% Generalized Aw-Rascal(GAR) model, with fundamental diagrams families has 3
% free parameters(TP), in data from California(CA)
% Riemann solver, HLL;
% fundamental diagram fitted from macroscopic data;
% initial density and boundary conditions are formed from NGSIM data
% Shimao Fan,
% Temple University, Math department,
% Fet 7th, 2012
%***********************************
% D === NGSIM_DensData
% V ==== NGSIM_VelData
% cl, cp and ca are parameters of fundamental diagrams
% h is size of grid point
% rhom is the maximum density
% n1 is average of time for interior data
% n2 is for boundary
% x0 is space vector choose h = 32 feet
% hom = 0 ..... homogeneous model
% hom = 1 ..... with stiff relaxiation terms
%***********************************
% Add relaxiation term on right hand side.
% Shimao Fan, March 15 2012
% Temple University.
%==========================================================================
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////
% Here FD is arctan function with three free parameters
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%------------------------------------------------
% determine dx and dt
% tfinal = 14;         % mins
h = x(2)-x(1);
N = length(x);       % length of x vector
k = 1/36000;         % length of time step, does not work start from ref = 8
t0 = 0:k:(tfinal/60);
%-------------------------------------------
dt = 2*k/ref; 
t = 0:dt:(tfinal/60);
D_lbc = spline(t0,BD(:,1),t);
D_rbc = spline(t0,BD(:,2),t);
V_lbc = spline(t0,BV(:,1),t);
V_rbc = spline(t0,BV(:,2),t);
% perform projection on the data
index_l = find(D_lbc>rhom);
index_r = find(D_rbc>rhom);
D_lbc(index_l) = .9*rhom;
V_lbc(index_l) = 1.0;
D_rbc(index_r) = .9*rhom;
V_rbc(index_r) = 1.0;
% interpol respect to time for all the data
for i = 1:size(D,2)   % loop respect space grids
    D_data(:,i) = spline(t0,D(:,i),t);
    V_data(:,i) = spline(t0,V(:,i),t);
end
%--------------------------------------------    
% at this end, we find time step satisfy CFL condition
M = ceil(tfinal/(60*dt));
lambda = dt/h;        % the length of road
%------------------------------------------------
% initial condition
icv = V(1,:);
icd = D(1,:);
index = find(icd>rhom);
icd(index) = .9*rhom;   % also do projection to initial data
icv(index) = 1.0;
U(1,:) = icd;
% calculate initial w vector by solving inverse problem
for i = 1:N
   U(2,i) = func_inverse_v_arctan_2p(icd(i),icv(i),q_min,q_midd,q_max,...
       ca,cb,rho0,v0,f0,uf,rhof,rhom,proj);
end

%*************************************************************************
%*******************  ----------------------  ****************************
%                       The main algorithm
%*******************  ----------------------  ****************************
%*************************************************************************

erho = zeros(1,M);
evel = zeros(1,M);
Q = unknownQ(U);
midd = round(N/2);   % middle position of study area

%--------------------------------

for n = 1:M      % loop over all time U = [density, w]; as unknowns

   %----------------------------------------------------------------------
   %  Insert new boundary conditions on right hand side
   rhol = D_lbc(n); rhor = D_rbc(n); % density of left and right bound
   vell = V_lbc(n); velr = V_rbc(n);
   %*****************************
   % apply bisection root finder
   %*****************************
   
   wl = func_inverse_v_arctan_2p(rhol,vell,q_min,q_midd,q_max,ca,cb,...
       rho0,v0,f0,uf,rhof,rhom,proj);
   wr = func_inverse_v_arctan_2p(rhor,velr,q_min,q_midd,q_max,ca,cb,...
       rho0,v0,f0,uf,rhof,rhom,proj);
   
   %**********************************
% %    % apply newton's iteration
   %**********************************
%     wl = w_max;
%     wr = wl;    % initial guess    
%     for i =1 : 2   % something wrong with this root finding process        
%        wl = wl+(vell-Velocity(rhol,wl,ca,cb,cc,rhom))./...
%            func_diff_v_w(rhol,wl,ca,cb,cc,rhom);
%        wr = wr+(velr-Velocity(rhor,wr,ca,cb,cc,rhom))./...
%           func_diff_v_w(rhor,wr,ca,cb,cc,rhom);       
%     end  
%     
%    wl = min(max(wl,w_min),w_max);
%    wr = min(max(wr,w_min),w_max); 
   %--------------------------------------
   lbc = [D_lbc(n);wl];
   rbc = [D_rbc(n);wr];
   Up1=[U(:,2:end),rbc];             % right boundary
   Um1=[lbc,U(:,1:end-1)];           % left boundary 

   Q = Q+lambda*(Num_Flux(Um1,U,ca,cb,rho0,v0,f0,uf,rhof,rhom)-...
       Num_Flux(U,Up1,ca,cb,rho0,v0,f0,uf,rhof,rhom));
   %----------------------------------------------------------------------
   %  find [rho,u] from conserved form, solution of next time step
   %-----------------------------------------------------------------
   U(1,:)=Q(1,:);                     % traffic density
   U(2,:)=Q(2,:)./Q(1,:);             % traffic velocity
   %----------------------------------------
   rho = U(1,:);
   q = U(2,:);
%    w = projection(w,w_min,w_max);
   vel = Velocity(rho,q,ca,cb,rho0,v0,f0,uf,rhof,rhom);
   
   num_rho(n) = rho(midd);
   num_vel(n) = vel(midd);
   %----------------------------------------
   
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
%    [e1,e2] = comput_error_nagim_time(D_data(n,:),...
%              V_data(n,:),rho,vel);
%    erho(n) = e1;
%    evel(n) = e2;
   e1 = norm(D_data(n,:)-rho,1);
   e2 = norm(V_data(n,:)-vel,1);
   erho(n) = e1;
   evel(n) = e2;
   % rescale with data
%    es1 = norm((D_data(n,:)-rho)./D_data(n,:),1);
%    es2 = norm((V_data(n,:)-vel)./V_data(n,:),1);
   
%    erhos(n) = es1;
%    evels(n) = es2;
%    n
end

% rescale of error
erho = erho/ref;
evel = evel/ref;


%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%---------------------------------------------------------
function y = func_diff_v_w(rho,w,ca,cb,cc,rhom)
     eps = 1.0e-06;
     dy = Velocity(rho,w+eps,ca,cb,cc,rhom) - ...
            Velocity(rho,w,ca,cb,cc,rhom);
     y = dy/eps;
%---------------------------------------------------------
% Ve function, the equalibrium function
function u = Ve(rho,lambda,p,alpha,rhom)  % density and 3 parameters   
      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = ((rho./rhom)-p).*lambda;
      f = alpha.*(a+(b-a).*(rho./rhom)-sqrt(1+y.^2));
      u = f./rho;    
%---------------------------------------------------------------
function  F = Num_Flux(Ul,Ur,ca,cb,rho0,v0,f0,uf,rhof,rhom) % numerical flux
   % find sending function
   s_flux = sending(Ul(1,:),Ul(2,:),ca,cb,rho0,v0,f0,uf,rhof,rhom);
   % receiving of flux
   r_flux = receiving(Ul,Ur,ca,cb,rho0,v0,f0,uf,rhof,rhom);
   % conservation of mass
   F(1,:) = min(s_flux,r_flux);   % take the minimum
   % conservation of momentum
   F(2,:) = Ul(2,:).*F(1,:);
%----------------------------------------------
% define the receiving function
function y = receiving(Ul,Ur,ca,cb,rho0,v0,f0,uf,rhof,rhom)
      w = Ul(2,:);
      
      a = polyval(ca,w); b = polyval(cb,w); 
%       a = ca(1)*w+ca(2);
%       b = cb(1)*w+cb(2);

      % determined by Ul
      rhoc = rho_critical_cgarz(a,b,rho0,v0,f0,rhom);  % a vector of rhoc
      % velocity on the right side
      ur = Velocity(Ur(1,:),Ur(2,:),ca,cb,rho0,v0,f0,uf,rhof,rhom);
      % newton iteration is faster/less costy
      rho_middle = Riemann_G(Ul,Ur,ca,cb,rho0,v0,f0,uf,rhof,rhom); % middle state
      
      case1 = rho_middle<=rhoc;  % maximum 
      case2 = rho_middle>rhoc;   % maximum potential
%       q_max = w;  % w represents the maximum flux
%       q_max = func_fd_seibold_2p(rhoc,a,b,rho0,v0,f0,uf,rhof,rhom); % maximum flow rate
      y = case2.*rho_middle.*ur + case1.*w;  
%----------------------------------------------
% define sending functions
function y = sending(rho,w,ca,cb,rho0,v0,f0,uf,rhof,rhom)% U = [\rho, w] two state variables

a = polyval(ca,w); b = polyval(cb,w); 
%       a = ca(1)*w+ca(2);
%       b = cb(1)*w+cb(2);
%       alpha = polyval(ca,w);     % alpha parameter
      rhoc = rho_critical_cgarz(a,b,rho0,v0,f0,rhom);  % a vector of rhoc/3-params fd
      case1 = rho<=rhoc;  % maximum possible
      case2 = rho>rhoc;   % maximum
%       q_max = w;   % w is the maximum flux
%       q_max = func_fd_seibold_2p(rhoc,a,b,rho0,v0,f0,uf,rhof,rhom); % maximum flow rate
      q = func_fd_seibold_2p(rho,a,b,rho0,v0,f0,uf,rhof,rhom);
      y = case1.*q + case2.*w;
%-----------------------------------------------
% calculate the middle density
function rho = Riemann_G(Ul,Ur,ca,cb,rho0,v0,f0,uf,rhof,rhom) % calculate the riemann solution/intermediate state
 % Ul and Ur are initial data of Riemann problem
  u = Velocity(Ur(1,:),Ur(2,:),ca,cb,rho0,v0,f0,uf,rhof,rhom); % right velocity
  w = Ul(2,:);     % left propoty
  rho = Ul(1,:);   % not sure which one is a good pick
  for i =1 : 2   % something wrong with this root finding process
       rho = rho+(u-Velocity(rho,w,ca,cb,rho0,v0,f0,uf,rhof,rhom))./...
           DiffVel(rho,w,ca,cb,rho0,v0,f0,uf,rhof,rhom);     
  end  
%--------------------------------------------------------
function y = unknownQ(U)    % the conserved unknowns [rho, rho*w]
        y(1,:) = U(1,:);
        y(2,:) = U(1,:).*U(2,:);
%----------------------------------------------------   
function u = DiffVel(rho,q,ca,cb,rho0,v0,f0,uf,rhof,rhom)% define derivative of velocity
       a = polyval(ca,q); b = polyval(cb,q);  % input data
%       a = ca(1)*q+ca(2);
%       b = cb(1)*q+cb(2);
      u = func_diff_vel_seibold_2p(rho,a,b,rho0,v0,f0,uf,rhof,rhom);
%----------------------------------------------------
function u = Velocity(rho,q,ca,cb,rho0,v0,f0,uf,rhof,rhom) % velocity function
       % rho is density;   
       a = polyval(ca,q); b = polyval(cb,q);
%       a = ca(1)*q+ca(2);
%       b = cb(1)*q+cb(2);
      u = func_fd_seibold_2p(rho,a,b,rho0,v0,f0,uf,rhof,rhom)./rho; 
%--------------------------------------------------------------
% function f = traffic_flux_smooth(rho,w,ca,cb,cc,rhom)% smooth flux function
%      % x_lambda and x_p are the coefficients of the fitting polynomials
%       a = polyval(ca,w); b = polyval(cb,w);  % input data
%       c = polyval(cc,w);
% %       a = cl(1)*w+cl(2);
% %       b = cp(1)*w+cp(2);
% %       c = ca(1)*w+ca(2);
%       f = func_fd_seibold_3p(rho,a,b,c,rhom);

%%%%%%%%%%%%%%%%%%%%%%End of Code%%%%%%%%%%%%%%%%%%%%%%

