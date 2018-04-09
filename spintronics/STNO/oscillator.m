clear;
%参数定义部分
alpha=0.005;     %α阻尼系数
gamma=2.2127614886*10^5;      %γ旋磁率
u0=4*pi*10^-7;           %真空磁导率
e=1.60217662*10^-19;      %电子电荷
hbar=1.05457266*10^-34;      %约化普朗克常数

H_initial=0.14/u0;     %有效场
H_ext=0;      %外加磁场

P=0.33;        %自旋极化率
Ms=13*10^5;    %饱和磁化率
t_lm=2*10^-9;     %自由层厚度
d_MgO=1.5*10^-9;    %MgO厚度
Area=pi*50*50*10^-18;   %有效面积
I=-3000*10^-6;        %电流
% J=0;
beta=-0.1;                 %field-like torque系数
VCMA=2*200*10^-15/(Ms*d_MgO*t_lm*u0);      %假设VCMA系数是200fJ/(A.m)

Hs_e=hbar*I*P/(2*e*u0*t_lm*Ms*Area);
step=1/1000;
N=2^24;
n=(0:N-1)/(N*step);

for i=1:1:2000/step
    time(i)=i*step;
    if i>50/step
        V(i)=2.5+0.5*cos(2*pi*0.2*step*i);
        %         V(i)=2+0.5;
    else
        V(i)=0;
    end
    %     if i==200/step
    %         I=-5000*10^-6;
    %         Hs_e=hbar*I*P/(2*e*u0*t_lm*Ms*Area);
    %     end
    if i==1
        H_eff(i)=H_initial;
        Hs(i)=Hs_e;
        theta(i)=0.1;
        phi(i)=0;
    else
        %                         Hs(i)=Hs_e*(2.25^2*P/(2.25^2+1+(2.25^2-1)*my(i-1)));
        H_eff(i)=H_initial+V(i)*VCMA;
        Hs(i)=Hs_e/(1+0.38*my(i-1));
        Hsf=Hs(i)*beta;
        %自己推导常用LLG方程结果
        theta(i)=theta(i-1)+(step*gamma*10^-9/(1+alpha^2))*(-H_eff(i)*alpha*cos(theta(i-1))*sin(theta(i-1))-...
            H_ext*alpha*sin(theta(i-1))+...
            Hsf*(-alpha*sin(phi(i-1))*cos(theta(i-1))-cos(phi(i-1)))+...
            Hs(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1))));
        phi(i)=phi(i-1)+(step*gamma*10^-9/(1+alpha^2))*(H_eff(i)*cos(theta(i-1))+H_ext+...
            Hsf*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1))+...
            Hs(i)*(cos(phi(i-1))+alpha*sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1)));
    end
    if theta(i)>pi
        theta(i)=2*pi-theta(i);
    end
    if theta(i)<0
        theta(i)=-theta(i);
    end
    mz(i)=cos(theta(i));
    mx(i)=sin(theta(i))*cos(phi(i));
    my(i)=sin(theta(i))*sin(phi(i));
%     if(rem(i,100)==0)
%         my_fft(round(i/100))=my(i);
%     end
end
figure;
subplot(2,1,1);
plot(time,theta);ylabel('theta');
subplot(2,1,2);
plot(time,phi);ylabel('phi');
figure;
subplot(3,1,1);
% plot(time,mx);ylabel('mx');
subplot(3,1,2);
% plot(time,mz);ylabel('mz');
plot(time(60/step:80/step),V(60/step:80/step));ylabel('V');
subplot(3,1,3);
plot(time(60/step:80/step),my(60/step:80/step));ylabel('my');
figure;
Y=fft(my(50/step:end),N);
% P=Y.*conj(Y)/N;
P=abs(Y);
P(1)=0;
plot(n(1:round(0.006*N)),P(1:round(0.006*N))); xlabel('Frequency');
max_p=max(P);
fre=find(P==max(P));
frequency=fre/(N*step)
% for j=1:length(theta)
%     x(j)=sin(theta(j))*cos(phi(j));
%     y(j)=sin(theta(j))*sin(phi(j));
%     z(j)=cos(theta(j));
% end
% figure;
% x1=load('C:\Users\Stuart\Documents\MATLAB\DTMF\x.txt');  %绘制一个标准的球形
% y1=load('C:\Users\Stuart\Documents\MATLAB\DTMF\y.txt');
% z1=load('C:\Users\Stuart\Documents\MATLAB\DTMF\z.txt');
% value1=spcrv([[x1(1) x1 x1(end)];[y1(1) y1 y1(end)];[z1(1) z1 z1(end)]],150);
% plot_std=plot3(value1(1,:),value1(2,:),value1(3,:));
% set(plot_std,'color',[180 180 180]/255);
% hold on;
%
% value1=spcrv([[x(1) x x(end)];[y(1) y y(end)];[z(1) z z(end)]],150);   %插值绘制翻转的曲线
% plot3(value1(1,:),value1(2,:),value1(3,:),'-k','linewidth',2);
