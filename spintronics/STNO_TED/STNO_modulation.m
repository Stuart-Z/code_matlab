%% 参数定义部分
clear;

alpha=0.005;     %α阻尼系数
gamma=2.2127614886*10^5;      %γ旋磁率
u0=4*pi*10^-7;           %真空磁导率
e=1.60217662*10^-19;      %电子电荷
hbar=1.05457266*10^-34;      %约化普朗克常数

H_k=2.09/u0;     %有效场        3mT时为2.094
H_ext=0.0/u0;      %外加磁场

P=0.50;        %自旋极化率
Ms=1.56*10^6;    %饱和磁化率
t_lm=2*10^-9;     %自由层厚度
d_MgO=1.5*10^-9;    %VCMA氧化层
Area=pi*150*50*10^-18/2;   %有效面积
I=-2800*10^-6;        %电流
beta=-0.1;
VCMA=2*200*10^-15/(Ms*d_MgO*t_lm*u0);      %假设VCMA系数是200fJ

if(I>-600*10^-6)               %模拟H_effective 随着电流会有一定的变化，根据Kubota论文 Kubota_2013_Appl._Phys._Express_6_103003
    H_eff_0=H_k*(1-0.082*(I/10^-3)^2)-Ms;
else
    H_eff_0=H_k*(1-0.082*(600*10^-6/10^-3)^2)-Ms;
end
%H_eff=H_k*(1-0.082*(I/10^-3)^2)-Ms;


% Ic=4*alpha*e*Ms*Area*t_lm*(H_eff+H_ext)/(hbar*(P^3));
% Ic

% Hs_e=hbar*I/(2*e*u0*t_lm*Ms*Area);
step=1/1000;        %仿真步长 单位ns
N=2^20;             %计算频率使用fft 定义fft为2^20
n=(0:N-1)/(N*step);
t_signal=10;    %定义码元长度  单位ns
t_stop=200;     %定义仿真跑多久

len=t_stop/step;
[I_sum, time, Hs, Hsf, theta, phi, mx, my, mz, V, Resistence, Voltage, ]=deal(zeros(1,len));
%% 电压配置
for i=1:1:t_stop/step   %生成VCMA电压的调制控制信号，在100ns后在开始调制
    V(i)=0;
    temp=(i-100/step)/(t_signal/step);
    temp1=floor(temp);
    temp2=mod(temp1,2);
    if temp2==0
        if((i*step-100-temp1*t_signal)<(t_signal/15))
            V(i)=0.8;        %0.91
        else
            V(i)=0.8;
        end
    else
        if((i*step-100-temp1*t_signal)<(t_signal/20))
            V(i)=1.25;        %1.14
        else
            V(i)=1.25;
        end
    end
end

%% LLG动力学仿真
for i=1:1:t_stop/step
    time(i)=i*step;
    I_sum(i)=I;
    Hs_e=hbar*I_sum(i)/(2*e*u0*t_lm*Ms*Area);
    
    H_eff=H_eff_0+V(i)*VCMA;
    if i==1
        Hs(i)=Hs_e;
        theta(i)=0.8;
        phi(i)=0;
    else
        Hs(i)=Hs_e*P/(1+P^2*my(i-1));
        Hsf(i)=Hs(i)*beta;
        %自己推导常用LLG方程结果
        theta(i)=theta(i-1)+(step*gamma*10^-9/(1+alpha^2))*(-H_eff*alpha*cos(theta(i-1))*sin(theta(i-1))-...
            H_ext*alpha*sin(theta(i-1))+...
            Hsf(i)*(-alpha*sin(phi(i-1))*cos(theta(i-1))-cos(phi(i-1)))+...
            Hs(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1))));
        phi(i)=phi(i-1)+(step*gamma*10^-9/(1+alpha^2))*(H_eff*cos(theta(i-1))+H_ext+...
            Hsf(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1))+...
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
    Resistence(i)=408+92*my(i);
    Voltage(i)=-I_sum(i)*Resistence(i);
end

%% fft计算频率
Y=fft(my(100/step:(100+t_signal)/step),N);
P=Y.*conj(Y)/N;
P(1)=0;
fre=find(P==max(P));
frequency1=fre/(N*step)

Y=fft(my((100+t_signal)/step:(100+2*t_signal)/step),N);
P=Y.*conj(Y)/N;
P(1)=0;
fre=find(P==max(P));
frequency2=fre/(N*step)

Y=fft(my((100+2*t_signal)/step:(100+3*t_signal)/step),N);
P=Y.*conj(Y)/N;
P(1)=0;
fre=find(P==max(P));
frequency3=fre/(N*step)

Y=fft(my(100/step:t_stop/step),N);
P=Y.*conj(Y)/N;
P(1)=0;
% figure;
% plot(n,P);

%% 画图显示效果
% figure;
% subplot(2,1,1);
% plot(time,theta);ylabel('theta');
% subplot(2,1,2);
% plot(time,phi);ylabel('phi');
% figure;
% subplot(3,1,1);
% plot(time(500/step:505/step),mx(500/step:505/step));ylabel('mx');
% subplot(3,1,2);
% plot(time(500/step:505/step),mz(500/step:505/step));ylabel('mz');
% subplot(3,1,3);
% plot(time(500/step:505/step),my(500/step:505/step));ylabel('my');
figure;
subplot(3,1,1);
plot(time,mx);ylabel('mx');
subplot(3,1,2);
plot(time,mz);ylabel('mz');
subplot(3,1,3);
plot(time,my);ylabel('my');
figure;
subplot(2,1,1)
plot(time,Voltage);ylabel('voltage');
subplot(2,1,2);
plot(time,V);ylabel('voltage dc');


