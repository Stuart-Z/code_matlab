alpha=0.01;
hs=0.1;
hp=0;
theta=[];
phi=[];
t_array=[];
mz=[];
mx=[];
x=[];
y=[];
z=[];
for i=1:800
    t_array(i)=i;
    if i==1
        theta(i)=0.01;
        phi(i)=pi/2;
    else
        theta(i)=theta(i-1)-sin(theta(i-1))*(alpha*cos(theta(i-1))+hp*(sin(phi(i-1))+alpha*cos(theta(i-1))*cos(phi(i-1)))*cos(phi(i-1))-hs);
        phi(i)=phi(i-1)-(cos(theta(i-1))+hp*(cos(phi(i-1))*cos(theta(i-1))-alpha*sin(phi(i-1)))*cos(phi(i-1))+hs*alpha);
    end
    if theta(i)>pi
        theta(i)=2*pi-theta(i);
    end
    if theta(i)<0
        theta(i)=-theta(i);
    end
    mz(i)=cos(theta(i));
    mx(i)=sin(theta(i))*cos(phi(i));
end
for j=1:600
    x(j)=sin(theta(j))*cos(phi(j));
    y(j)=sin(theta(j))*sin(phi(j));
    z(j)=cos(theta(j));
end
figure;
plot(t_array,theta,'-');
figure;
plot(t_array,mz,'-');
figure;
[x1,y1,z1]=sphere;
% mesh(x1,y1,z1);
% hold on;
plot3(x,y,z,'-');