function SNR_test( SNR )
alpha=0.005;     %������ϵ��
gamma=2.2127614886*10^5;      %��������
u0=4*pi*10^-7;           %��մŵ���
e=1.60217662*10^-19;      %���ӵ��
hbar=1.05457266*10^-34;      %Լ�����ʿ˳���

H_k=2.09/u0;     %��Ч��        3mTʱΪ2.094
H_ext=0.0/u0;      %��Ӵų�

P=0.50;        %����������
Ms=1.56*10^6;    %���ʹŻ���
t_lm=2*10^-9;     %���ɲ���
d_MgO=1.5*10^-9;
Area=pi*150*50*10^-18/2;   %��Ч���
I=-1600*10^-6;        %����
beta=-0.1;
V=0;
VCMA=2*200*10^-15/(Ms*d_MgO*t_lm*u0);      %����VCMAϵ����200fJ

if(I>-600*10^-6)
    H_eff_0=H_k*(1-0.082*(I/10^-3)^2)-Ms;
else
    H_eff_0=H_k*(1-0.082*(600*10^-6/10^-3)^2)-Ms;
end
% H_eff=H_k*(1-0.082*(I/10^-3)^2)-Ms;


% Ic=4*alpha*e*Ms*Area*t_lm*(H_eff+H_ext)/(hbar*(P^3));
% Ic

% Hs_e=hbar*I/(2*e*u0*t_lm*Ms*Area);
step=1/1000;
N=2^20;
n=(0:N-1)/(N*step);
t_signal=125;    %������Ԫ����  ��λns
f_high=2;
f_low=1.6;
% SNR=63;

%% �����ź�ģ��
% for i=1:1:1500/step
%     time(i)=i*step;
%     N(i)=SNR;
%     f_ac_GHz(i)=1;
%     if i>100/step
%         if i>100/step&&i<=500/step
%             f_ac_GHz(i)=2;
%         end
%         if i>500/step&&i<=1000/step
%             f_ac_GHz(i)=f_low;
%         end
%         if i>1000/step&&i<1500/step
%             f_ac_GHz(i)=2;
%         end
%         I_ac(i)=150*10^-6*cos(2*pi*step*i*f_ac_GHz(i));   %RF�źŵķ���Ϊx00uA
%     else
%         I_ac(i)=0;
%     end
% end
for i=1:1:600/step
    time(i)=i*step;
    N(i)=SNR;
     f_ac_GHz(i)=f_low;
    temp=(i-100/step)/(t_signal/step);
    temp1=floor(temp);
    temp2=mod(temp1,2);
    if(i>100/step)
        if temp2==0
            f_ac_GHz(i)=f_high;
        else
            f_ac_GHz(i)=f_low;
        end
        I_ac(i)=150*10^-6*cos(2*pi*step*i*f_ac_GHz(i));   %RF�źŵķ���Ϊx00uA
    else
        I_ac(i)=0;
    end
end

%% �����ͨ�˲���
% figure
% plot(time,I_ac);
I_ac=awgn(I_ac,SNR);        %���������

len=length(time);
X=fft(I_ac,len);           %�����ͨ�˲���ʵ�֣�����fft���Ƶ�� �����ȫ��Ϊ��
% figure;
% plot(time,X);title('�˲�ǰƵ��')
for i=1:1:len
    if(i>ceil(len/490)&&i<ceil(len*(1-1/490)))   %����1/2T�Ľ�ֹƵ����һ����ͨ��
        X(i)=0;
    end
    if(i<ceil(len*f_low/(510*2))||i>ceil(len*(1-f_low/(510*2))))   %����1/2T�Ľ�ֹƵ����һ����ͨ�˲�   ��Ϊ1Gʱ��1/1100     1.5GʱΪ1/700  1.8GΪ560 
        X(i)=0;
    end    
end
% figure;
% plot(time,X);title('�˲���Ƶ��');
magX=abs(X);
angX=angle(X);
Y=magX.*exp(1i*angX);
I_ac=ifft(Y,len);         %������Ҷ�任���µõ��ź� ��ʱ�ź��൱���˲����
% figure;
% plot(time,I_ac);title('�˲����ź�');

%% LLG���̶���ѧģ��
for i=1:1:600/step
    time(i)=i*step;
    I_sum(i)=I_ac(i)+I;
    Hs_e=hbar*I_sum(i)/(2*e*u0*t_lm*Ms*Area);
    H_eff=H_eff_0+V*VCMA;
    if i==100/step
        V=0.5;
    end
%     if i==1100/step
%         V=1.03;
%     end
%     if i==2100/step
%         V=0.48;
%     end
    if i==1
        Hs(i)=Hs_e;
        theta(i)=0.8;
        phi(i)=0;
    else
        Hs(i)=Hs_e*P/(1+P^2*my(i-1));
        Hsf(i)=Hs(i)*beta;
        %�Լ��Ƶ�����LLG���̽��
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
% figure;
% plot(time,I_sum);title('�˲����ź�');

%% �����ͨ
len=length(time);
X=fft(Voltage,len);           %�����ͨ�˲���ʵ�֣�����fft���Ƶ�� �����ȫ��Ϊ��
% figure;
% plot(time,X);
for i=1:1:len
    if(i>ceil(len/10000)&&i<ceil(len*9999/10000))   %����1/2T�Ľ�ֹƵ����һ����ͨ�˲�
        X(i)=0;
    end
end
% figure;
% plot(time,X);
magX=abs(X);
angX=angle(X);
Y=magX.*exp(1i*angX);
Voltage_dc=ifft(Y,len);         %������Ҷ�任���µõ��ź� ��ʱ�ź��൱���˲����

%% ��ͼ��ʾЧ��
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
figure;
subplot(3,1,1)
plot(time,Voltage);ylabel('voltage');
subplot(3,1,2);
plot(time(200/step:end),Voltage_dc((200/step:end)));ylabel('voltage dc');
subplot(3,1,3);
plot(time,N);



end

