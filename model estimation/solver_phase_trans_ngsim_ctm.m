function [erho,evel] = solver_phase_trans_ngsim_ctm(tfinal,D,V,BD,BV,...
    x,ref,u_r,rho_r,q_min,q_max,rhom)
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
% close all;
% some parameters
% pick the middle curve/ NGSIM data
% lda = P(1);
% p = P(2);
% alpha = P(3);
% rhovec = .1:.1:rhom;           % traffic density vector       
%------------------------------------------------
% determine dx and dt
% tfinal = 14;         % mins
h = x(2)-x(1);
N = length(x);       % length of x vector
k = 1/36000;         % length of time step, does not work start from ref = 8
t0 = 0:k:(tfinal/60);
m1 = length(t0);
% sl = 80;
%-------------------------------------------
dt = 2*k/ref; 
t = 0:dt:(tfinal/60);
D_lbc = spline(t0,BD(1:m1,1),t);
D_rbc = spline(t0,BD(1:m1,2),t);
V_lbc = spline(t0,BV(1:m1,1),t);
V_rbc = spline(t0,BV(1:m1,2),t);
% interpol respect to time for all the data
for i = 1:size(D,2)   % loop respect space grids
     D_data(:,i) = spline(t0,D(1:m1,i),t);
     V_data(:,i) = spline(t0,V(1:m1,i),t);
end
%--------------------------------------------    
% at this end, we find time step satisfy CFL condition
M = ceil(tfinal/(60*dt));
lambda = dt/h;        % the length of road
% L = x(end);      
% pos = ceil(.224/h);
%------------------------------------------------
% initial condition
icv = V(1,:);
icd = D(1,:);
U(1,:) = icd;
% icw = 0*icv;
% find out w data
% for n = 1:length(icv)
%    rho = icd(n); vel = icv(n);       % density of left and right bound
%    drho = abs(rho-rhovec); 
%    [C,indrho] = min(drho);
%    uvec = squeeze(Mw(indrho,2,:));   % velocity vector   
%    dvel = abs(vel-uvec);
%    wvec = squeeze(Mw(indrho,1,:));   % w vector
%    [C,indv] = min(dvel); 
%    icw(n) = wvec(indv); 
% end
% whos
% icbeta = 0*icv+1;
% 
% for i =1 : 2   % something wrong with this root finding process       
%   icbeta = icbeta+(icv-Velocity(icd,icbeta,ct,sl,rhom))./...
%            func_diff_v_beta(icd,icbeta,ct,sl,rhom);
% end

for i = 1:length(icd)
  U(2,i) = func_inverse_v_pase_trans(icd(i),icv(i),q_min,q_max,u_r,rho_r,rhom);
end

%*************************************************************************
%*******************  ----------------------  ****************************
%                       The main algorithm
%*******************  ----------------------  ****************************
%*************************************************************************

erho = zeros(1,M);
evel = zeros(1,M);
% w_min = 55;
% w_max = 80;
Q = unknownQ(U);     % unknowns of in conserved form [rho, flux]
midd = round(N/2);   % middle position of study area

for n = 1:M      % loop over all time U = [density, w]; as unknowns
   % represent Q as function of U, the conservative unknowns
%    Q(1,:)=U(1,:);
%    Q(2,:)=Q(1,:).*U(2,:);
   %----------------------------------------------------------------------
   %  Insert new boundary conditions on right hand side
   rhol = D_lbc(n); rhor = D_rbc(n); % density of left and right bound
   ul = V_lbc(n); ur = V_rbc(n);
   
   % apply bisection root finder
   ql = func_inverse_v_pase_trans(rhol,ul,q_min,q_max,u_r,rho_r,rhom);
   qr = func_inverse_v_pase_trans(rhor,ur,q_min,q_max,u_r,rho_r,rhom);


   lbc = [D_lbc(n);ql];
   rbc = [D_rbc(n);qr];
   Up1=[U(:,2:end),rbc];             % right boundary
   Um1=[lbc,U(:,1:end-1)];           % left boundary 

   %---------------------------------------------------
   Q=Q+lambda*(Num_Flux(Um1,U,u_r,rho_r,rhom)- ...
       Num_Flux(U,Up1,u_r,rho_r,rhom));   
   %----------------------------------------------------------------------
   %  find [rho,u] from conserved form, solution of next time step
   %-----------------------------------------------------------------
   U(1,:)=Q(1,:);                     % traffic density
   U(2,:)=Q(2,:)./Q(1,:);             % traffic velocity
   %----------------------------------------
   rho = U(1,:);
   q = U(2,:);
%    beta = min(max(beta,beta_min),beta_max);

%    w = projection(w,w_min,w_max);
   vel = func_u_phase_transition(rho,q,u_r,rho_r,rhom);
   
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
   e1 = norm(D_data(n,:)-rho,1);
   e2 = norm(V_data(n,:)-vel,1);
   erho(n) = e1;
   evel(n) = e2;
end

%---------------------------------
% erho = erho/ref/rhom;
% evel = evel/ref/vm;

erho = erho/ref;
evel = evel/ref;

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
        
function  F = Num_Flux(Ul,Ur,u_r,rho_r,rhom) % numerical flux
   % find sending function
   s_flux = sending(Ul(1,:),Ul(2,:),u_r,rho_r,rhom);
   % receiving of flux
   r_flux = receiving(Ul,Ur,u_r,rho_r,rhom);
   % conservation of mass
   F(1,:) = min(s_flux,r_flux);   % take the minimum
   % conservation of momentum
   F(2,:) = Ul(2,:).*F(1,:);
%----------------------------------------------
% define the receiving function
function y = receiving(Ul,Ur,u_r,rho_r,rhom)
      w = Ul(2,:); % determined by Ul
      rhoc = rho_critical(w,u_r,rho_r);  % a vector of rhoc
      % velocity on the right side
      ur = func_u_phase_transition(Ur(1,:),Ur(2,:),u_r,rho_r,rhom);
      % newton iteration is faster/less costy
      rho_middle = Riemann_G(Ul,Ur,u_r,rho_r,rhom); % middle state
      
      case1 = rho_middle<=rhoc;  % maximum 
      case2 = rho_middle>rhoc;   % maximum potential
      q_max = traffic_flux(rhoc,w,u_r,rho_r,rhom); % maximum flow rate
%       q_max = func_fd_seibold_2p(rhoc,a,b,rho0,v0,f0,uf,rhof,rhom); % maximum flow rate
      y = case2.*rho_middle.*ur + case1.*q_max;  
%----------------------------------------------
% define sending functions
function y = sending(rho,w,u_r,rho_r,rhom)% U = [\rho, w] two state variables
      rhoc = rho_critical(w,u_r,rho_r);  % a vector of rhoc/3-params fd
      case1 = rho<=rhoc;  % maximum possible
      case2 = rho>rhoc;   % maximum
%       q_max = w;   % w is the maximum flux
      q_max = traffic_flux(rhoc,w,u_r,rho_r,rhom); % maximum flow rate
      q = traffic_flux(rho,w,u_r,rho_r,rhom);
      y = case1.*q + case2.*q_max;
%-----------------------------------------------
% calculate the intermediate state density
function rho = Riemann_G(Ul,Ur,u_r,rho_r,rhom) % calculate the riemann solution/intermediate state
 % Ul and Ur are initial data of Riemann problem
  u = func_u_phase_transition(Ur(1,:),Ur(2,:),u_r,rho_r,rhom); % right velocity
  w = Ul(2,:);     % left propoty
  rho = Ul(1,:);   % not sure which one is a good pick
  for i =1 : 2   % something wrong with this root finding process
       rho = rho+(u-func_u_phase_transition(rho,w,u_r,rho_r,rhom))./...
           DiffVel(rho,w,u_r,rho_r,rhom);     
  end   
%---------------------------------------------
function y = unknownQ(U)    % the conserved unknowns [rho, rho*beta]
        y(1,:) = U(1,:);
        y(2,:) = U(1,:).*U(2,:);
%-----------------------------------------------------
function y = rho_critical(q,u_r,rho_r)
    y = rho_r/2.*(1-sqrt(1-4*q./(u_r.*rho_r)));
%-----------------------------------------------------        
function y = traffic_flux(rho,q,u_r,rho_r,rhom)
     u = func_u_phase_transition(rho,q,u_r,rho_r,rhom);
     y = rho.*u;
%----------------------------------------------------   
function y = DiffVel(rho,q,u_r,rho_r,rhom)% define derivative of velocity
      eps = 1.0e-05;
      y =(func_u_phase_transition(rho+eps,q,u_r,rho_r,rhom)-...
          func_u_phase_transition(rho,q,u_r,rho_r,rhom))/eps;
%---------------------------------------------------------------- 
function y = minimum(a,b) % here a and b could be vectors
    temp1 = a<=b;
    temp2 = a>b;
    y = temp1.*a+temp2.*b;    % where y is the minimum
  
function y = maximum(a,b) % so as above function
    temp1 = a<=b;
    temp2 = a>b;
    y = temp1.*b+temp2.*a;

%%%%%%%%%%%%%%%%%%%%%%End of Code%%%%%%%%%%

