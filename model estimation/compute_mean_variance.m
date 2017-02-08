clear;
mean_er=zeros(8,1);
mean_ev=zeros(8,1);
var_er=zeros(8,1);
var_ev=zeros(8,1);


i_end=74;


load data_dens_minn_2.mat             % sensor 2 data from 4pm tp 5pm/D
D = Dens;
D_mean=mean(D');
D_index=zeros(i_end,1);

for i=1:i_end
    if D_mean(1,i)>80
        D_index(i,1)=1;
    end
end

count0=0;
for i =1:i_end  
    if rem(i-1,6)~=0 && D_index(i,1)==1
        count0=count0+1;
    end
end
    
N=count0;
count_day= zeros(N,1);
count1=0;
for i =1:i_end  
    if rem(i-1,6)~=0 && D_index(i,1)==1
        count1=count1+1;
        count_day(count1,1)=i;
    end
end


e_rho=zeros(N,8);
e_vel=zeros(N,8);
count=0;
for i =1:i_end
    
    if rem(i-1,6)~=0 && D_index(i,1)==1
        count=count+1;
        
        [erho,evel]=comput_minn_validation_cluster_newptfd(i);
        e_rho(count,:)=erho;
        e_vel(count,:)=evel;
%         for j=1:8
%             mean_er(j)=mean_er(j)+(1/N)*erho(1,j);
%             mean_ev(j)=mean_ev(j)+(1/N)*evel(1,j);
%         end
    end
end

% for j=1:8
%     for i =1:N
%             var_er(j)=var_er(j)+(1/(N-1))*(e_rho(i,j)-mean_er(j))*(e_rho(i,j)-mean_er(j));
%             var_ev(j)=var_ev(j)+(1/(N-1))*(e_vel(i,j)-mean_ev(j))*(e_vel(i,j)-mean_ev(j));
%         
%     end
% end
