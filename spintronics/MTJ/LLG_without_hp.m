global alpha;
alpha=0.01;
global hs;
hs=-0.02;
n_array=[];
mz=[];
mx=[];
x=[];
y=[];
z=[];
for i=1:490
    n_array(i)=i;
    mz(i)=cos(f1(i));
    mx(i)=sin(f1(i))*cos(f2(i));
end
for j=1:360
    x(j)=sin(f1(j))*cos(f2(j));
    y(j)=sin(f1(j))*sin(f2(j));
    z(j)=cos(f1(j));
end
figure;
plot(n_array,mz,'-');
figure;
plot(n_array,mx,'-');
figure;
plot3(x,y,z,'-');