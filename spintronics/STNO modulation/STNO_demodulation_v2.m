%% �������岿��
clear;
alpha=0.005;     %������ϵ��
gamma=2.2127614886*10^5;      %��������
u0=4*pi*10^-7;           %��մŵ���
e=1.60217662*10^-19;      %���ӵ��
hbar=1.05457266*10^-34;      %Լ�����ʿ˳���

H_k=2.09/u0;     %��Ч��        3mTʱΪ2.094
H_ext=0.0/u0;      %��Ӵų�

P_factor=0.50;        %����������
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
% n=(0:N-1)/(N*step);
t_signal=100;                                                       %������Ԫ����  ��λns
t_stop=100+t_signal*280;                                            %�������ĳ���
signal_begin=40*t_signal+100;
signal_end=signal_begin+201*t_signal;
f_high=2;
f_low=1.6;
SNR=42;

%% �����ź�ģ��
len=t_stop/step;
[f_ac_GHz, I_ac, time]=deal(zeros(1,len));
for i=1:1:t_stop/step
    time(i)=i*step;
    temp=(i-100/step)/(t_signal/step);
    temp1=floor(temp);
    temp2=mod(temp1,2);
    if temp2==0
        f_ac_GHz(i)=f_high;
    else
        f_ac_GHz(i)=f_low;
    end
    I_ac(i)=-180*10^-6*cos(2*pi*step*i*f_ac_GHz(i));   %RF�źŵķ���Ϊx00uA
end
clear f_ac_GHz;
%% �����ͨ�˲���
I_ac=awgn(I_ac,SNR);        %���������
% figure
% plot(time,I_ac);

X=fft(I_ac,len);           %�����ͨ�˲���ʵ�֣�����fft���Ƶ�� �����ȫ��Ϊ��
% figure;
% Power=X.*conj(X)/N;
% Power(1)=0;
% plot(time,Power);title('�˲�ǰƵ��')
Y=zeros(1,len);
for i=1:1:len
    if((i>ceil(len*(f_low-1/t_signal)/(1000))&&i<ceil(len*(f_high+1/t_signal)/1000))||(i>ceil(len*(1-(f_high+1/t_signal)/1000))&&i<ceil(len*(1-(f_low-1/t_signal)/(1000)))))
        Y(i)=X(i);
    end
end
% figure;
% Power=Y.*conj(Y)/N;
% Power(1)=0;
% plot(time,Power);title('�˲���Ƶ��');
magX=abs(Y);
angX=angle(Y);
Y=magX.*exp(1i*angX);
I_ac=ifft(Y,len);         %������Ҷ�任���µõ��ź� ��ʱ�ź��൱���˲����
I_ac=real(I_ac);
% figure;
% plot(time,I_ac);title('�˲����ź�');
clear X Y angX magX;

%% LLG���̶���ѧģ�� STNO_1 2GHz
[Hs, Hsf, I_sum, mx, my, mz, phi, Resistence, theta, Voltage1]=deal(zeros(1,len));
for i=1:1:t_stop/step
    time(i)=i*step;
    I_sum(i)=I_ac(i)+I;
    Hs_e=hbar*I_sum(i)/(2*e*u0*t_lm*Ms*Area);
    H_eff=H_eff_0+V*VCMA;
    if i>100/step
        V=0.51;
    end
    if i==1
        Hs(i)=Hs_e;
        theta(i)=0.8;
        phi(i)=0;
    else
        Hs(i)=Hs_e*P_factor/(1+P_factor^2*my(i-1));
        Hsf(i)=Hs(i)*beta;
%         �Լ��Ƶ�����LLG���̽��
        theta(i)=theta(i-1)+(step*gamma*(10^-9)/(1+alpha^2))*(-H_eff*alpha*cos(theta(i-1))*sin(theta(i-1))-...
            H_ext*alpha*sin(theta(i-1))+...
            Hsf(i)*(-alpha*sin(phi(i-1))*cos(theta(i-1))-cos(phi(i-1)))+...
            Hs(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1))));
        phi(i)=phi(i-1)+(step*gamma*(10^-9)/(1+alpha^2))*(H_eff*cos(theta(i-1))+H_ext+...
            Hsf(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1))+...
            Hs(i)*(cos(phi(i-1))+alpha*sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1)));
    end
    mz(i)=cos(theta(i));
    mx(i)=sin(theta(i))*cos(phi(i));
    my(i)=sin(theta(i))*sin(phi(i));
    Resistence(i)=408+92*my(i);
    Voltage1(i)=-I_sum(i)*Resistence(i);
end
clear Hs Hsf I_sum phi Resistence theta ;
clear mx my mz;
% %% fft����Ƶ��
% Y=fft(my((signal_begin)/step:(signal_begin+t_signal)/step),N);
% P=Y.*conj(Y)/N;
% P(1)=0;
% fre=find(P==max(P));
% frequency1=fre/(N*step)
% 
% Y=fft(my((signal_begin+t_signal)/step:(signal_begin+2*t_signal)/step),N);
% P=Y.*conj(Y)/N;
% P(1)=0;
% fre=find(P==max(P));
% frequency2=fre/(N*step)
% 
% Y=fft(my((signal_begin+2*t_signal)/step:(signal_begin+3*t_signal)/step),N);
% P=Y.*conj(Y)/N;
% P(1)=0;
% fre=find(P==max(P));
% frequency3=fre/(N*step)
%% �����ͨ
len=length(time);
X=fft(Voltage1,len);           %�����ͨ�˲���ʵ�֣�����fft���Ƶ�� �����ȫ��Ϊ��
% figure;
% Power=X.*conj(X)/N;
% Power(1)=0;
% plot(time,Power);
Y=zeros(1,len);
for i=1:1:len
    if(i<ceil(len/(1000*t_signal))||i>ceil(len*(1-1/(1000*t_signal))))   %����1/2T�Ľ�ֹƵ����һ����ͨ�˲�
        Y(i)=X(i);
    end
end
Y(1)=0;Y(2)=0;
% figure;
% Power=Y.*conj(Y)/N;
% Power(1)=0;
% plot(time,Power);
magX=abs(Y);
angX=angle(Y);
Y=magX.*exp(1i*angX);
Voltage_dc1=ifft(Y,len);         %������Ҷ�任���µõ��ź� ��ʱ�ź��൱���˲����
Voltage_dc1=real(Voltage_dc1);
clear X Y angX magX;

%% LLG���̶���ѧģ�� STNO_2 1.6GHz
[Hs, Hsf, I_sum, mx, my, mz, phi, Resistence, theta, Voltage2]=deal(zeros(1,len));
for i=1:1:t_stop/step
    I_sum(i)=I_ac(i)+I;
    Hs_e=hbar*I_sum(i)/(2*e*u0*t_lm*Ms*Area);
     H_eff=H_eff_0+V*VCMA;
    if i>100/step
        V=0.31;
    end
    if i==1
        Hs(i)=Hs_e;
        theta(i)=0.8;
        phi(i)=0;
    else
        Hs(i)=Hs_e*P_factor/(1+P_factor^2*my(i-1));
        Hsf(i)=Hs(i)*beta;
        %�Լ��Ƶ�����LLG���̽��
        theta(i)=theta(i-1)+(step*gamma*(10^-9)/(1+alpha^2))*(-H_eff*alpha*cos(theta(i-1))*sin(theta(i-1))-...
            H_ext*alpha*sin(theta(i-1))+...
            Hsf(i)*(-alpha*sin(phi(i-1))*cos(theta(i-1))-cos(phi(i-1)))+...
            Hs(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1))));
        phi(i)=phi(i-1)+(step*gamma*(10^-9)/(1+alpha^2))*(H_eff*cos(theta(i-1))+H_ext+...
            Hsf(i)*(-alpha*cos(phi(i-1))+sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1))+...
            Hs(i)*(cos(phi(i-1))+alpha*sin(phi(i-1))*cos(theta(i-1)))/sin(theta(i-1)));
    end
    mz(i)=cos(theta(i));
    mx(i)=sin(theta(i))*cos(phi(i));
    my(i)=sin(theta(i))*sin(phi(i));
    Resistence(i)=408+92*my(i);
    Voltage2(i)=-I_sum(i)*Resistence(i);
end
clear Hs Hsf I_sum phi Resistence theta ;
clear mx my mz;
% %% fft����Ƶ��
% Y=fft(my((signal_begin)/step:(signal_begin+t_signal)/step),N);
% P=Y.*conj(Y)/N;
% P(1)=0;
% fre=find(P==max(P));
% frequency4=fre/(N*step)
% 
% Y=fft(my((signal_begin+t_signal)/step:(signal_begin+2*t_signal)/step),N);
% P=Y.*conj(Y)/N;
% P(1)=0;
% fre=find(P==max(P));
% frequency5=fre/(N*step)
% 
% Y=fft(my((signal_begin+2*t_signal)/step:(signal_begin+3*t_signal)/step),N);
% P=Y.*conj(Y)/N;
% P(1)=0;
% fre=find(P==max(P));
% frequency6=fre/(N*step)
%% �����ͨ
len=length(time);
X=fft(Voltage2,len);           %�����ͨ�˲���ʵ�֣�����fft���Ƶ�� �����ȫ��Ϊ��
% figure;
% Power=X.*conj(X)/N;
% Power(1)=0;
% plot(time,Power);
Y=zeros(1,len);
for i=1:1:len
    if(i<ceil(len/(1000*t_signal))||i>ceil(len*(1-1/(1000*t_signal))))   %����1/2T�Ľ�ֹƵ����һ����ͨ�˲�
        Y(i)=X(i);
    end
end
Y(1)=0;Y(2)=0;
% figure;
% Power=Y.*conj(Y)/N;
% Power(1)=0;
% plot(time,Power);
magX=abs(Y);
angX=angle(Y);
Y=magX.*exp(1i*angX);
Voltage_dc2=ifft(Y,len);         %������Ҷ�任���µõ��ź� ��ʱ�ź��൱���˲����
Voltage_dc2=real(Voltage_dc2);
clear X Y angX magX;


%% �����о�ͼ��
[th, voltage_recover, signal_input]=deal(zeros(1,len));

threshold=mean(Voltage_dc1((signal_begin)/step:signal_end/step));      %������ֵ��ѹ
% th(1)=threshold;                                                    %Ϊ����ͼ����ʾ��ֵ��ѹ ���Ǽ�����һ������
signal_input(1)=0;                                                  %����Ķ������ź�
voltage_recover(1)=0;                                               %�ָ��Ķ������ź�
number_error=0;                                                     %��������Ĵ���  �ܴ���
number_total=0;
number_error_0to1=0;
number_error_1to0=0;
for i=2:1:t_stop/step
    voltage_recover(i)=voltage_recover(i-1);                        %ֻ���о���仯
%     th(i)=threshold;
    if(mod((i-signal_begin/step-t_signal/(2*step)),t_signal/step)==0&&i>(signal_begin/step)&&i<(signal_end/step))
        if(Voltage_dc1(i)>Voltage_dc2(i))                                 %����ֵ�Ƚ�
            voltage_recover(i)=0;
        else
            voltage_recover(i)=1;
        end
    end
    temp=(i-100/step-t_signal/(2*step))/(t_signal/step);            %�����ź�����ʱ�ļ���
    temp1=floor(temp);
    temp2=mod(temp1,2);
    if temp2==0
        signal_input(i)=1;
    else
        signal_input(i)=0;
    end
    if(mod((i-signal_begin/step),t_signal/step)==0&&i>(signal_begin/step)&&i<(signal_end/step))
        if(voltage_recover(i)~=signal_input(i))                     %�жϻָ�����ȷ��
            number_error=number_error+1;
            if(voltage_recover(i)==1)
                number_error_0to1=number_error_0to1+1;
            else
                number_error_1to0=number_error_1to0+1;
            end
        end
        number_total=number_total+1;
    end
end

%% ��ͼ��ʾЧ��
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

% figure;
% subplot(3,1,1);
% plot(time,mx);ylabel('mx');
% subplot(3,1,2);
% plot(time,mz);ylabel('mz');
% subplot(3,1,3);
% plot(time,my);ylabel('my');

figure;
% subplot(2,1,1)
% % plot(time,Voltage1);ylabel('voltage');
% subplot(2,1,2);
% plot(time(200/step:end),Voltage_dc1(200/step:end));ylabel('voltage dc');
% hold on;
% plot(time(200/step:end),Voltage_dc2(200/step:end),'r');
% figure;
% plot(time(200/step:end),1.2*voltage_recover(200/step:end));ylabel('voltage recover');
% hold on;
% plot(time(200/step:end),signal_input(200/step:end),'r');ylabel('voltage input');


