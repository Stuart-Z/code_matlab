%% Writen by Zuodong Zhang (2017/12/26)%%%
clear;
tic;
global gamma  alpha H_VCMA_dv Hext Hk Ms...
     kB T mu0 Hs_dI mp

%% LLG parameters %%
%%%------------
q=1.6e-19;                              % Coulombs
hbar=6.626e-34/2/pi;                    % Reduced Planck's constant (J-s)
alpha = 0.1;                            %0.01 Gilbert damping parameter
g = 1.76e7*4*pi/1000;                   %1.76e7 Gyromagnetic ratio [(rad)/(Oe.s)]
mu0=1.2566e-6;                          %Magnetic permeability
uB=9.27400949e-24;
T=300;
kB=1.38e-23;

%% Magnet Parameters (taken from experiment)

Ms = 1.2e5; % A/m, Saturation Magnetization [emu/cm^3] emu/cm^3 = 1 KA/m
Ku2 = 6624; % Uni. anisotropy constant [J/m^3]
Hk = 2*Ku2/mu0/Ms; % Switching field [A/m]
Hk;
%%%%%%%%%%%%%     VCMA[MKS]    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%72.864%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=0.33;
V_vcma=1.0598;                        %1.05984 HVCMA=HK Voltage across the MTJ,[V]
ksi=10e-15;                    %VCMA constant [J/V.m];37e-10[erg/V.cm]
t_fl=1e-9;                      %thickness of freelayer
Area=16*16*10^-18;
Vol=t_fl*Area;                  % Volume [m^3]
gamma=g/(1+alpha^2);            % reduced factor
d_MgO=1.6e-9;                   % thickness of MgO barrier[m]
Hs_dI=hbar*P/(2*q*mu0*t_fl*Ms*Area);
H_VCMA_dv=2*ksi/(mu0*Ms*d_MgO*t_fl); %*mz,[Oe]
H_VCMA=V_vcma*H_VCMA_dv;
H_VCMA;
delta=(mu0*(Hk-H_VCMA)*Ms*Vol)/(2*kB*T);
delta;
% I_c=(alpha*g*q/uB/P)*Ms*Vol*(Hk-H_VCMA);    %from zhao Recent progresses in STT-MRAM 2016
I_c=2*q*alpha*2*delta*kB*T/hbar/P;        %from datta 2017
I=2*10^-6;                            %current

Hext=0;                         % external field
mp=[0 0 1];                     % pinned layer
theta0=sqrt(kB*T/mu0/Ms/Hk/Vol);  % initial angle
step=1e-2;

%% Initial conditions of the simulation %%
mz = cos(theta0); % Magnet slightly off easy axis
m_initial = [0 sqrt(1-mz^2) mz];

%%   LLG     %%
time1 = 5000; % duration ns
tot1 = time1/step;
delta_t=step*10^-9; % time step
[theta,phi,Time,Energy]=deal(zeros(1,tot1));
[m1,m2,m3]=deal(zeros(tot1,3));
m1_positive=0;  m2_positive=0;  m3_positive=0;
m1_negtive=0;   m2_negtive=0;   m3_negtive=0;
choice=1;
for i=1:1:tot1
    Time(i)=i*step;
 
    %discription of thermal noise
    for j=1:1:3
        r=0+1.*randn(1,3); % Normal distribution, expectation = 0, standard deviation = 1
        Hth(j,:)=(sqrt(2*kB*T*alpha/(mu0*Ms*gamma*Vol*delta_t)))*r; % thermal noise field
    end

    if i==1
        [m1(i,:),m2(i,:),m3(i,:),]=deal(m_initial);
    else

%         I1=I*(-2+2*m2_read+2*m3_read);
%         I2=I*(1+2*m1_read-m3_read);
%         I3=I*(1+2*m1_read-m2_read);
        if(rem(i,250/step)==0)
            choice=~choice;
        end
        
        if (choice==1)
            I1=3*I*(m2_read);
        else
            I1=3*I*(m3_read);
        end
        I2=1*10^-6;
        I3=-2*10^-6;
        
        m1(i,:)=sLLG_step(I1, m1(i-1,:), Hth(1,:), V_vcma, delta_t);
        m2(i,:)=sLLG_step(I2, m2(i-1,:), Hth(2,:), V_vcma, delta_t);
        m3(i,:)=sLLG_step(I3, m3(i-1,:), Hth(3,:), V_vcma, delta_t);

    end
    m1_read=sign(m1(i,3));    m2_read=sign(m2(i,3));    m3_read=sign(m3(i,3));   
    

    %count
    if(m1(i,3)>0)
        m1_positive=m1_positive+1;
    else
        m1_negtive=m1_negtive+1;
    end
    
    if(m2(i,3)>0)
        m2_positive=m2_positive+1;
    else
        m2_negtive=m2_negtive+1;
    end
    if(m3(i,3)>0)
        m3_positive=m3_positive+1;
    else
        m3_negtive=m3_negtive+1;
    end
end

%% plot
probability1=m1_positive/(m1_positive+m1_negtive);
probability2=m2_positive/(m2_positive+m2_negtive);
probability3=m3_positive/(m3_positive+m3_negtive);
figure;
subplot(3,1,1);
plot(Time,m1(:,3),'-','linewidth',1.0); % m_z
title('m1');
subplot(3,1,2);
plot(Time,m2(:,3),'-','linewidth',1.0); % m_z
title('m2');
subplot(3,1,3);
plot(Time,m3(:,3),'-','linewidth',1.0); % m_z
title('m3');
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%