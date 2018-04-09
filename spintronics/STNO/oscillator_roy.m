function  frequency  = oscillator_roy( I )
%参数定义部分
alpha=0.01;     %α阻尼系数
gamma=2.2127614886*10^5;      %γ旋磁率
u0=4*pi*10^-7;           %真空磁导率
e=1.60217662*10^-19;      %电子电荷
hbar=1.05457266*10^-34;      %约化普朗克常数

H_k=1.4266/u0;     %有效场
H_ext=0/u0;      %外加磁场

H_eff=H_k-1.3823/u0;
% H_eff=H_k*(1-0.012*(I/300)^2)-1.3823/u0;

P=0.57;        %自旋极化率
Ms=11*10^5;    %饱和磁化率
t_lm=1.6*10^-9;     %自由层厚度
Area=pi*35*75*10^-18;   %有效面积
% if I<500
%     beta=-0.25*((620-I/10)*I/(570*500))^2;
% else
%     beta=-0.25;
% end
I=-I*10^-6;        %电流
beta=-0.2;

% Ic=4*alpha*e*Ms*Area*t_lm*(H_eff+H_ext)/(hbar*(P^3));
% Ic

Hs_e=hbar*I/(2*e*u0*t_lm*Ms*Area);
step=1/1000;
N=2^20;
n=(0:N-1)/(N*step);

for i=1:1:1500/step
    time(i)=i*step;
    if i==1
        Hs(i)=Hs_e;
        theta(i)=1;
        phi(i)=0;
    else
%           Hs(i)=Hs_e*(2.25^2*P/(2.25^2+1+(2.25^2-1)*my(i-1)));
%         Hs(i)=Hs_e*P/(1+0.325*my(i-1));
        Hs(i)=Hs_e;
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
    
end
% figure;
% subplot(2,1,1);
% plot(time,theta);ylabel('theta');
% subplot(2,1,2);
% plot(time,phi);ylabel('phi');
% figure;
% subplot(3,1,1);
% plot(time,mx);ylabel('mx');
% subplot(3,1,2);
% plot(time,mz);ylabel('mz');
% subplot(3,1,3);
% plot(time,my);ylabel('my');
% figure;
Y=fft(my,N);
P=Y.*conj(Y)/N;
P(1)=0;
% plot(n(1:8000),P(1:8000));
max_p=max(P);
fre=find(P==max(P));
frequency=fre/(N*step)
end

