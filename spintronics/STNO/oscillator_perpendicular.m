function  frequency  = oscillator_perpendicular( I , V)
%�������岿��
alpha=0.005;     %������ϵ��
gamma=2.2127614886*10^5;      %��������
u0=4*pi*10^-7;           %��մŵ���
e=1.60217662*10^-19;      %���ӵ��
hbar=1.05457266*10^-34;      %Լ�����ʿ˳���

H_eff=0.14/u0;     %��Ч��
H_ext=0/u0;      %��Ӵų�

P=0.33;        %����������
Ms=13*10^5;    %���ʹŻ���
t_lm=2*10^-9;     %���ɲ���
d_MgO=1.5*10^-9;
Area=pi*50*50*10^-18;   %��Ч���
I=-I*10^-6;        %����
% J=0;
beta=-0.1;
VCMA=2*200*10^-15/(Ms*d_MgO*t_lm*u0);      %����VCMAϵ����200fJ

H_eff=H_eff*(1-0.18*(I/1000)^2);
Hs_e=hbar*I*P/(2*e*u0*t_lm*Ms*Area);
step=1/1000;
N=2^20;
n=(0:N-1)/(N*step);

for i=1:1:1500/step
    time(i)=i*step;
    if i==100/step
        H_eff=H_eff+V*VCMA;
    end
    if i==1
        Hs(i)=Hs_e;
        theta(i)=0.2;
        phi(i)=pi;
    else
%                         Hs(i)=Hs_e*(2.25^2*P/(2.25^2+1+(2.25^2-1)*my(i-1)));
        Hs(i)=Hs_e/(1+0.38*my(i-1));
        Hsf=Hs(i)*beta;
        %����SUN��ģ��
%         theta(i)=theta(i-1)+(step*gamma*10^-9/(1+alpha^2))*(-H_eff*alpha*cos(theta(i-1))*sin(theta(i-1))-...
%             H_ext*alpha*sin(theta(i-1))+...
%             Hsf*(alpha*sin(phi(i-1))*cos(theta(i-1))-cos(phi(i-1)))+...
%             Hs(i)*(alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1))));
%         phi(i)=phi(i-1)+(step*gamma*10^-9/(1+alpha^2))*(-H_eff*cos(theta(i-1))-H_ext+...
%             Hsf*(alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1))+...
%             Hs(i)*(cos(phi(i-1))-alpha*sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1)));
        %�Լ��Ƶ�����LLG���̽��
        theta(i)=theta(i-1)+(step*gamma*10^-9/(1+alpha^2))*(-H_eff*alpha*cos(theta(i-1))*sin(theta(i-1))-...
            H_ext*alpha*sin(theta(i-1))+...
            Hsf*(-alpha*sin(phi(i-1))*cos(theta(i-1))-cos(phi(i-1)))+...
            Hs(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1))));
        phi(i)=phi(i-1)+(step*gamma*10^-9/(1+alpha^2))*(H_eff*cos(theta(i-1))+H_ext+...
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
% figure;
Y=fft(my,N);
P=Y.*conj(Y)/N;
P(1)=0;
% plot(n(1:20000),P(1:20000));
max_p=max(P);
fre=find(P==max(P));
frequency=fre/(N*step)
end

