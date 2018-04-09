alpha=0.1;
hs=0.76;
hp=4.8;
theta=[];
phi=[];
t_array=[];
t_array2=[];
mz=[];
mx=[];
x=[];
y=[];
z=[];
for i=1:100
    t_array2(i)=i;
    if i==1
        theta(i)=0.01;
        phi(i)=pi/2;
    else
        for k=0:100
            theta(i+k)=theta(i+k-1)-0.1*(sin(theta(i+k-1))*(alpha*cos(theta(i+k-1))+hp*(sin(phi(i+k-1))+alpha*cos(theta(i+k-1))*cos(phi(i+k-1)))*cos(phi(i+k-1))-hs));
            phi(i+k)=phi(i+k-1)-0.1*((cos(theta(i+k-1))+hp*(cos(phi(i+k-1))*cos(theta(i+k-1))-alpha*sin(phi(i+k-1)))*cos(phi(i+k-1))+hs*alpha));
            t_array(i+k)=i+k;
        end
        theta(i)=theta(i+10);
        phi(i)=phi(i+10);
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
% for j=1:150
%     x(j)=sin(theta(j))*cos(phi(j));
%     y(j)=sin(theta(j))*sin(phi(j));
%     z(j)=cos(theta(j));
% end
figure;
plot(t_array2,mz,'-');
figure;
plot(t_array2,mx,'-');
% figure;
% % [x1,y1,z1]=sphere;
% % mesh(x1,y1,z1);
% % hold on;
% plot3(x,y,z,'-');