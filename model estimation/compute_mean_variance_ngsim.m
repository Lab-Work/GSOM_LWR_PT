clear;
error_m=zeros(6,7);
k=3;


count=0;
for i =1:3
    i
    for j=6:6
        j
        [erho,evel]=comput_ngsim_validation_cgarz_paper(i,j,k);
        error_m((i-1)*2+1,j)=erho;
        error_m((i-1)*2+2,j)=evel;
    end
end

