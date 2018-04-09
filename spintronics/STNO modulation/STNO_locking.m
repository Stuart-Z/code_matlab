function [frequency,fre_predict,mz_mean,mz_cal] = STNO_locking(I,V,f_ac_GHz )
%%%% STNO器件特性仿真函数 输入VCMA电压、直流电流、交流小信号频率  
%%%%　输出仿真的频率的mz均值以及根据energy balance theory计算的频率和mz
%%参数定义部分
alpha=0.005;     %α阻尼系数
gamma=2.2127614886*10^5;      %γ旋磁率
u0=4*pi*10^-7;           %真空磁导率
e=1.60217662*10^-19;      %电子电荷
hbar=1.05457266*10^-34;      %约化普朗克常数

H_k=2.09/u0;     %有效场       
H_ext=0.0/u0;      %外加磁场

eta=0.50;        %自旋极化率
Ms=1.56*10^6;    %饱和磁化率
t_lm=2*10^-9;     %自由层厚度
d_MgO=1.5*10^-9;    %VCMA的氧化层厚度
Area=pi*150*50*10^-18/2;   %有效面积
I=-I*10^-6;        %电流
beta=-0.1;          %field like torque coefficient
% V=0;
VCMA=2*200*10^-15/(Ms*d_MgO*t_lm*u0);      %假设VCMA系数是200fJ

% H_eff=H_k-Ms;
if(I>-600*10^-6)            %模拟H_effective 随着电流会有一定的变化，根据Kubota论文 Kubota_2013_Appl._Phys._Express_6_103003
    H_eff=H_k*(1-0.082*(I/10^-3)^2)-Ms;
else
    H_eff=H_k*(1-0.082*(600*10^-6/10^-3)^2)-Ms;
end

% Hs_e=hbar*I/(2*e*u0*t_lm*Ms*Area);
step=1/1000;        %仿真步长 单位ns
N=2^20;             %计算频率使用fft 定义fft为2^20
n=(0:N-1)/(N*step);
len=1500/step;
[I_sum, time, Hs, Hsf, theta, phi, mx, my, mz ]=deal(zeros(1,len));
%% 交流信号生成
for i=1:1:1500/step
    if(i<300/step)
        I_sum(i)=I;
    else
        I_sum(i)=I+00*10^-6*cos(2*pi*step*f_ac_GHz*i);   %加入的交流信号   幅度300uA
    end
%     I_ac(i)=100*10^-6*cos(2*pi*step*f_ac_GHz*i);
end

%% 利用LLG方程对器件动力过程仿真
for i=1:1:1500/step
    time(i)=i*step;
    Hs_e=hbar*I_sum(i)/(2*e*u0*t_lm*Ms*Area);  %电流自旋转移矩大小
    if i==100/step                             %在100ns后等起振且基本稳定再加入VCMA电压
        H_eff=H_eff+V*VCMA;
    end
    if i==1
        Hs(i)=Hs_e;
        theta(i)=0.2;
        phi(i)=0;
    else
        Hs(i)=Hs_e*eta/(1+eta^2*my(i-1));
        Hsf(i)=Hs(i)*beta;                  %field like torque
        %自己推导常用LLG方程结果 my方向pinned layer
        theta(i)=theta(i-1)+(step*gamma*10^-9/(1+alpha^2))*(-H_eff*alpha*cos(theta(i-1))*sin(theta(i-1))-...
            H_ext*alpha*sin(theta(i-1))+...
            Hsf(i)*(-alpha*sin(phi(i-1))*cos(theta(i-1))-cos(phi(i-1)))+...
            Hs(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1))));
        phi(i)=phi(i-1)+(step*gamma*10^-9/(1+alpha^2))*(H_eff*cos(theta(i-1))+H_ext+...
            Hsf(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1))+...
            Hs(i)*(cos(phi(i-1))+alpha*sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1)));
    end
    mz(i)=cos(theta(i));
    mx(i)=sin(theta(i))*cos(phi(i));
    my(i)=sin(theta(i))*sin(phi(i));
    mz_cal(i)=(1-(Hs(i)/(alpha*H_eff))^2)^0.5;
end
%% 能量平衡理论 公式来自Theoretical Study of Spin-Torque Oscillator with Perpendicularly Magnetized Free Layer
for i=0.001:0.001:1.57
    theta_cal(round(i*1000))=i;
%     current(round(i*1000))=2*alpha*e*u0*eta*Ms*t_lm*Area*(H_ext+H_eff*cos(i))*sin(i)*sin(i)/(hbar*cos(i)*(1/(1-eta^4*sin(i)*sin(i))^0.5-1));
%     error_i(round(i*1000))=abs(current(round(1000*i))+I);
    W_s(round(i*1000))=hbar*eta*abs(I)*(H_ext+H_eff*cos(i))*cos(i)*(1/(1-eta^4*sin(i)*sin(i))^0.5-1)/(2*e*eta^2*Ms*t_lm*Area*u0)+...
        beta*(hbar*eta*abs(I)/(2*e*eta^2*Ms*t_lm*Area*u0))^2*((1+eta^4*cos(2*i))/(1-eta^4*sin(i)*sin(i))^1.5-1);
    W_a(round(i*1000))=alpha*(H_ext+H_eff*cos(i))^2*sin(i)*sin(i)+...
        alpha*(hbar*beta*eta*abs(I)/(2*e*eta^2*Ms*t_lm*Area*u0))^2*((1+eta^4*cos(2*i))/(1-eta^4*sin(i)*sin(i))^1.5-1)+...
        2*alpha*beta*(hbar*eta*abs(I)/(2*e*eta^2*Ms*t_lm*Area*u0))*(H_ext+H_eff*cos(i))*cos(i)*(1/(1-eta^4*sin(i)*sin(i))^0.5-1);
    error_i(round(i*1000))=abs(W_a(round(i*1000))-W_s(round(i*1000)));
end
% figure;
% plot(theta_cal,current);
% hold on;
% plot(theta_cal,error_i,'r');
mz_cal=cos(find(error_i==min(error_i(800:end)))/1000);  %寻找两个能量相平衡的点
% figure;
% plot(theta_cal,W_s);
% hold on;
% plot(theta_cal,W_a,'r');

%% FFT计算频率
% figure;
Y=fft(my(300/step:end),N);
% Y=fft(my,N);
P=Y.*conj(Y)/N;
Power=sum(P);
P(1)=0;
% plot(n,P);
fre=find(P==max(P));
frequency=fre/(N*step);
frequency;
% Power
mz_mean=mean(mz(400/step:1000/step));

% fre_predict=gamma*(H_ext+(H_eff)*(mz_mean-0.015))/(2*pi*10^9);
fre_predict=gamma*(H_ext+(H_eff)*(mz_cal))/(2*pi*10^9);
fre_predict

% figure;
% subplot(3,1,1);
% plot(time,mx);ylabel('mx');
% hold on;
% % % plot(time,8*10^3*I_ac,'r');
% subplot(3,1,2);
% plot(time,mz);ylabel('mz');
% hold on;
% plot(time,mz_cal,'r');
% subplot(3,1,3);
% plot(time,my);ylabel('my');
% hold on;
% plot(time,4*10^3*I_ac,'r');

end

