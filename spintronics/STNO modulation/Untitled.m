clear;


%参数定义部分
alpha=0.005;     %α阻尼系数
gamma=2.2127614886*10^5;      %γ旋磁率
u0=4*pi*10^-7;           %真空磁导率
e=1.60217662*10^-19;      %电子电荷
hbar=1.05457266*10^-34;      %约化普朗克常数

H_k=2.09/u0;     %有效场       
H_ext=0.00/u0;      %外加磁场

eta=0.50;        %自旋极化率
Ms=1.56*10^6;    %饱和磁化率
t_lm=2*10^-9;     %自由层厚度
d_MgO=1.5*10^-9;
Area=pi*150*50*10^-18/2;   %有效面积
I=-800*10^-6;        %电流
beta=-0.1;
V=0;
VCMA=2*200*10^-15/(Ms*d_MgO*t_lm*u0);      %假设VCMA系数是200fJ

if(I>-600*10^-6)
    H_eff=H_k*(1-0.082*(I/10^-3)^2)-Ms;
else
    H_eff=H_k*(1-0.082*(600*10^-6/10^-3)^2)-Ms;
end

% H_eff=H_k-Ms;
% Hs_e=hbar*I/(2*e*u0*t_lm*Ms*Area);
step=1/1000;
N=2^20;
n=(0:N-1)/(N*step);
len=1500/step;
[I_sum, time, Hs, Hsf, theta, phi, mx, my, mz ]=deal(zeros(1,len));
for i=1:1:1500/step
    if(i<300/step)
        I_sum(i)=I;
    else
        I_sum(i)=I;
    end
%     I_ac(i)=100*10^-6*cos(2*pi*step*f_ac_GHz*i);
end


for i=1:1:1500/step
    time(i)=i*step;
    Hs_e=hbar*I_sum(i)/(2*e*u0*t_lm*Ms*Area);
    if i==100/step
        H_eff=H_eff+V*VCMA;
    end
    if i==1
        Hs(i)=Hs_e;
        theta(i)=0.2;
        phi(i)=0;
    else
        Hs(i)=Hs_e*eta/(1+eta^2*my(i-1));
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
    mz(i)=cos(theta(i));
    mx(i)=sin(theta(i))*cos(phi(i));
    my(i)=sin(theta(i))*sin(phi(i));
end
%% 能量平衡理论
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
mz_cal1=cos(find(error_i==min(error_i(800:end)))/1000);
figure;
plot(theta_cal,W_s);
hold on;
plot(theta_cal,W_a,'r');

for i=1:1:1500/step
    mz_cal(i)=mz_cal1;
end

%% FFT计算频率
% figure;
Y=fft(my(300/step:end),N);
% Y=fft(my,N);
P=Y.*conj(Y)/N;
% P(1)=0;
Power=sum(P);
P(1)=0;
% plot(n,P);
fre=find(P==max(P));
frequency=fre/(N*step);
frequency;
% Power
mz_mean=mean(mz(400/step:1000/step));
% figure;
% plot(time,Hs);  
% figure;
% subplot(2,1,1);
% plot(time,theta);ylabel('theta');
% subplot(2,1,2);
% plot(time,phi);ylabel('phi');

% fre_predict=gamma*(H_ext+(H_eff)*(mz_mean-0.015))/(2*pi*10^9);
fre_predict=gamma*(H_ext+(H_eff)*(mz_cal1))/(2*pi*10^9);
fre_predict
w0=gamma*(H_eff+H_ext);
sigma=gamma*hbar*eta^3/(4*e*u0*Ms*t_lm*Area);
I_threshold=alpha*w0/sigma;
I_threshold

figure;
subplot(3,1,1);
plot(time,mx);ylabel('mx');
% plot(time,8*10^3*I_ac,'r');
subplot(3,1,2);
plot(time,mz);ylabel('mz');
hold on;
plot(time,mz_cal,'r');
subplot(3,1,3);
plot(time,my);ylabel('my');
% plot(time,4*10^3*I_ac,'r');