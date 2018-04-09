alpha=0.01;
hs=-0.7;
hp=0;
h_ext=0;
theta=[];
phi=[];
t_array=[];
t_array2=[];
mz=[];
mx=[];
x=[];
y=[];
z=[];
for i=1:10000
    t_array2(i)=i;
    if i==1
        theta(i)=0.4;
        phi(i)=pi/2;
    else
        for k=0:100
            theta(i+k)=theta(i+k-1)-0.01*(sin(theta(i+k-1))*(alpha*cos(theta(i+k-1))+hp*(sin(phi(i+k-1))+alpha*cos(theta(i+k-1))*cos(phi(i+k-1)))*cos(phi(i+k-1))*sin(theta(i+k-1))+hs*sin(theta(i+k-1))+h_ext*alpha*sin(theta(i+k-1))));
            phi(i+k)=phi(i+k-1)-0.01*((cos(theta(i+k-1))+hp*(cos(phi(i+k-1))*cos(theta(i+k-1))-alpha*sin(phi(i+k-1)))*cos(phi(i+k-1))-hs*alpha+h_ext));
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
for j=1:10000
    x(j)=sin(theta(j))*cos(phi(j));
    y(j)=sin(theta(j))*sin(phi(j));
    z(j)=cos(theta(j));
end
% for j=1:1000
%     y(j)=sin(theta(j))*cos(phi(j));
%     z(j)=sin(theta(j))*sin(phi(j));
%     x(j)=cos(theta(j));
% end
figure;
value=spcrv([[t_array2(1) t_array2 t_array2(end)];[mz(1) mz mz(end)]],100);
%plot(t_array2,mz,'-');
plot(value(1,:)/10,value(2,:),'-');
title('magnetization switching');
xlabel('nomalization time');
ylabel('M_z');
% figure;
% plot(t_array2,mx,'-');
figure;
% value1=spcrv([[x(1) x x(end)];[y(1) y y(end)];[z(1) z z(end)]],150);
% plot3(value1(1,:),value1(2,:),value1(3,:),'-k','linewidth',2);
% hold on;
x1=load('C:\Users\Stuart\Documents\MATLAB\DTMF\x.txt');
y1=load('C:\Users\Stuart\Documents\MATLAB\DTMF\y.txt');
z1=load('C:\Users\Stuart\Documents\MATLAB\DTMF\z.txt');
value1=spcrv([[x1(1) x1 x1(end)];[y1(1) y1 y1(end)];[z1(1) z1 z1(end)]],150);
plot_std=plot3(value1(1,:),value1(2,:),value1(3,:));
set(plot_std,'color',[180 180 180]/255);
hold on;
value1=spcrv([[x(1) x x(end)];[y(1) y y(end)];[z(1) z z(end)]],150);
plot3(value1(1,:),value1(2,:),value1(3,:),'-k','linewidth',2);
% dlmwrite('C:\Users\Stuart\Documents\MATLAB\DTMF\x.txt',x);
% dlmwrite('C:\Users\Stuart\Documents\MATLAB\DTMF\y.txt',y);
% dlmwrite('C:\Users\Stuart\Documents\MATLAB\DTMF\z.txt',z);