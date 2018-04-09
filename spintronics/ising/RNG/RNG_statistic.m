clear all;
clc
tic


A=1000;
mz=[];
cont=0;


for i=1:A
    mz(i,1)=statistic(i);
end

for j=1:A
    if mz(j,1)<-0.95
        cont=cont+1;
    else cont=cont;
    end
end

p=cont/A;
% pp(jj,1)=p;

% end
toc