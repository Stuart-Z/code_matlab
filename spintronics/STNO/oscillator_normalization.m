clear;
alpha=0.01;
hs=0.9;
h=1;
step=1/100;
N=2^20;
n=(0:N-1)*2*pi/(N);

for i=1:1:400000
    time(i)=i*step;
    if i==1
        theta(i)=0.1;
        phi(i)=0;
    else
        theta(i)=theta(i-1)+step*(-alpha*sin(theta(i-1))*cos(theta(i-1))-h*alpha*sin(theta(i-1))+...
            hs*alpha*cos(phi(i-1))+hs*sin(phi(i-1))*cos(theta(i-1)));
        phi(i)=phi(i-1)+step*(-cos(theta(i-1))-h+hs*(cos(phi(i-1))-alpha*sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1)));
        %         theta(i)=theta(i-1)+0.001*(-alpha*sin(theta(i-1))*cos(theta(i-1))-h*alpha*sin(theta(i-1))-hs*sin(theta(i-1)));
        %         phi(i)=phi(i-1)+0.001*(-cos(theta(i-1))-h+hs*alpha);
    end
    if theta(i)>pi
        theta(i)=2*pi-theta(i);
    end
    if theta(i)<0
        theta(i)=-theta(i);
    end
    %         if i>1
    %             hs=hs-(theta(i)-theta(i-1));
    %         end
    mz(i)=cos(theta(i));
    mx(i)=sin(theta(i))*cos(phi(i));
    my(i)=sin(theta(i))*sin(phi(i));
    if mod(i,100)==0
        mx_fft(i/100)=mx(i);
    end
    
end
figure;
subplot(2,1,1);
plot(time,theta);ylabel('theta');
subplot(2,1,2);
plot(time,phi);ylabel('phi');
figure;
subplot(3,1,1);
plot(time,mx);ylabel('mx');
subplot(3,1,2);
plot(time,mz);ylabel('mz');
subplot(3,1,3);
plot(time,my);ylabel('my');
figure;
Y=fft(mx_fft,N);
P=Y.*conj(Y)/N;
plot(n,P);
for j=1:length(theta)
    x(j)=sin(theta(j))*cos(phi(j));
    y(j)=sin(theta(j))*sin(phi(j));
    z(j)=cos(theta(j));
end
figure;
x1=load('C:\Users\Stuart\Documents\MATLAB\DTMF\x.txt');  %绘制一个标准的球形
y1=load('C:\Users\Stuart\Documents\MATLAB\DTMF\y.txt');
z1=load('C:\Users\Stuart\Documents\MATLAB\DTMF\z.txt');
value1=spcrv([[x1(1) x1 x1(end)];[y1(1) y1 y1(end)];[z1(1) z1 z1(end)]],150);
plot_std=plot3(value1(1,:),value1(2,:),value1(3,:));
set(plot_std,'color',[180 180 180]/255);
hold on;

value1=spcrv([[x(1) x x(end)];[y(1) y y(end)];[z(1) z z(end)]],150);   %插值绘制翻转的曲线
plot3(value1(1,:),value1(2,:),value1(3,:),'-k','linewidth',2);
hold on;
value1=spcrv([[x(1) x x(end)];[y(1) y y(end)];[z(1) z z(end)]],150);
plot3(value1(1,:),value1(2,:),value1(3,:),'-k','linewidth',2);
