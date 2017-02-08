% function [num_rho,num_vel,erho,evel] = solver_arzq_ngsim(tfinal,D,V,BD,BV,x,ref,vm0,vm,rhom)
function [num_rho,num_vel,erho,evel] = solver_arzq_ngsim_ctm(tfinal,D,V,BD,BV,x,ref,vm,rhom)

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
% vm0 calculate the error, act as a weight, which is const for all model
% vm is to study the effect of model on different maximum velocity
%***********************************
% Add relaxiation term on right hand side.
% Shimao Fan, March 15 2012
% Temple University.
%==========================================================================
% close all;
% some parameters
% pick the middle curve/ NGSIM data
%------------------------------------------------
% determine dx and dt
% tfinal = 14;         % mins
h = x(2)-x(1);
N = length(x);       % length of x vector
midd = round(N/2);
k = 1/36000;         % length of time step, does not work start from ref = 8
t0 = 0:k:(tfinal/60);
m1 = length(t0);
%-------------------------------------------------
% dt = 2*k/ref; 
dt = k/ref; 

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
% if ref < 8
%    inter = 2^2/ref;
%    D_data = D(1:inter:end,:);
%    V_data = V(1:inter:end,:);
% else
%     D_data = D;
%     V_data = V;
% end
%--------------------------------------------    
% at this end, we find time step satisfy CFL condition
M = ceil(tfinal/(60*dt));
lambda = dt/h;        % the length of road
% L = x(end);      
% pos = ceil(.224/h);
%------------------------------------------------
% initial conditions
U(1,:) = D(1,:);
U(2,:) = V(1,:);              % initial condition
%*************************************************************************
%*******************  ----------------------  ****************************
%                       The main algorithm
%*******************  ----------------------  ****************************
%*************************************************************************

erho = zeros(1,M);
evel = zeros(1,M);

Q = unknownQ(U,vm,rhom);
   
for n = 1:M      % loop over all time U = [rho, v]; as unknowns
   % represent Q as function of U, the conservative unknowns
%    Q = unknownQ(U,vm,rhom);
   %----------------------------------------------------------------------
   %  Insert new boundary conditions on right hand side
   rhol = D_lbc(n); rhor = D_rbc(n); % density of left and right bound
   vell = V_lbc(n); velr = V_rbc(n);
   lbc = [rhol;vell];
   rbc = [rhor;velr];
   Up1=[U(:,2:end),rbc];    % right boundary
   Um1=[lbc,U(:,1:end-1)];  % left boundary 
   Q = Q+lambda*(Num_Flux(Um1,U,vm,rhom)-Num_Flux(U,Up1,vm,rhom));
   %----------------------------------------------------------------------
   %  find [rho,u] from conserved form, solution of next time step
   %-----------------------------------------------------------------
   U(1,:)=Q(1,:);                     % traffic density
   U(2,:)=Q(2,:)./Q(1,:)+Ve(U(1,:),vm,rhom)-vm;             % traffic velocity
   %----------------------------------------
   rho = U(1,:);
   vel = U(2,:);
   num_rho(n) = rho(midd);
   num_vel(n) = vel(midd);
%    U(:,1) = lbc;
%    U(:,end) = rbc;
   %----------------------
   %   The plots part
%    ----------------------
%   if mod(n,10) == 0
%    subplot(2,1,1)
%    plot(x,U(1,:),'m-','linewidth',2.5), hold on
%    plot(x,D(n,1:N),'b-','linewidth',2.5), hold off
%    title(sprintf('The traffic density at time t=%4.3f',n*dt))
%    axis([0,x(end),0,600])
%    ylabel('Density','fontsize',14)
%    xlabel('Space','fontsize',14)
%    set(gca,'linewidth',2)
%    set(gca,'fontsize',14)
% %    %----------------------------------------------------------------
%    subplot(2,1,2)
%    plot(x,U(2,:),'m-','linewidth',2.5),hold on
%    plot(x,V(n,1:N),'b-','linewidth',2.5), hold off
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
   
%       % rescale with data
%    es1 = norm((D_data(n,:)-rho)./D_data(n,:),1);
%    es2 = norm((V_data(n,:)-vel)./V_data(n,:),1);
%    
%    erhos(n) = es1;
%    evels(n) = es2;
end

erho = erho/ref;
evel = evel/ref;

% erhos = erhos/ref;
% evels = evels/ref;


%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
function  F = Num_Flux(Ul,Ur,vm,rhom) % numerical flux
   % find sending function
   w = Ul(2,:)+vm-Ve(Ul(1,:),vm,rhom);
   s_flux = sending(Ul(1,:),w,vm,rhom);
   % receiving of flux
   r_flux = receiving(Ul,Ur,vm,rhom);
   % conservation of mass
   F(1,:) = min(s_flux,r_flux);   % take the minimum
   % conservation of momentum
   F(2,:) = w.*F(1,:);
%----------------------------------------------
% define the receiving function
function y = receiving(Ul,Ur,vm,rhom)
      w = Ul(2,:)+vm-Ve(Ul(1,:),vm,rhom);
      % determined by Ul
      rhoc = rho_critical_arz(w,vm,rhom);  % a vector of rhoc/3-params fd
      % velocity on the right side
      ur = Ur(2,:);
      % newton iteration is faster/less costy
      rho_middle = Riemann_G(Ul,Ur,vm,rhom); % middle state
      
      case1 = rho_middle<=rhoc;  % maximum 
      case2 = rho_middle>rhoc;   % maximum potential
      
      q_max = traffic_flux(rhoc,w,vm,rhom); % maximum flow rate
      y = case2.*rho_middle.*ur + case1.*q_max;  
%----------------------------------------------
% define sending functions
function y = sending(rho,w,vm,rhom)% U = [\rho, w] two state variables
      rhoc = rho_critical_arz(w,vm,rhom);  % a vector of rhoc/3-params fd
      case1 = rho<=rhoc;  % maximum possible
      case2 = rho>rhoc;   % maximum
      q_max = traffic_flux(rhoc,w,vm,rhom); % maximum flow rate
      q = traffic_flux(rho,w,vm,rhom);
      y = case1.*q + case2.*q_max;
%-----------------------------------------------
% calculate the middle density
function rho = Riemann_G(Ul,Ur,vm,rhom) % calculate the riemann solution/intermediate state
 % Ul and Ur are initial data of Riemann problem
   % calculat w value
  ur = Ur(2,:);
  ul = Ul(2,:);
  ve = Ve(Ul(1,:),vm,rhom); % right velocity
  rho = rhom/vm*(vm+ul-ur-ve);
%---------------------------------------------
function y = unknownQ(U,vm,rhom)    % the conserved unknowns [rho, rho*(v-Ve)]
        y(1,:) = U(1,:);     % density
        y(2,:) = U(1,:).*(U(2,:)+vm-Ve(U(1,:),vm,rhom));
%----------------------------------------------------         
function rhoc = rho_critical_arz(w,vm,rhom)% equilibrium velocity  
     rhoc = w*rhom/(2*vm);      
%--------------------------------------------------------  
 function y = traffic_flux(rho,w,vm,rhom)% equilibrium velocity  
     u = Velocity(rho,w,vm,rhom);
     y = u.*rho;
%---------------------------------------------------------     
 function y = Velocity(rho,w,vm,rhom)% equilibrium velocity  
      ve = Ve(rho,vm,rhom); 
      y = (w-vm) + ve;
%--------------------------------------------------------   
function y = diffvel(vm,rhom)% define derivative of velocity
        y = -vm/rhom;
%----------------------------------------------------------------
% modify here
function y = Ve(rho,vm,rhom)  % density and 3 parameters 
     y = vm*(1-(rho/rhom));   % which is the quadratic form
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


