%% Writen by Zuodong Zhang (2017/12/09)%%%
clear;
tic;
global gamma  alpha H_VCMA_dv Hext Hk Ms...
     kB T mu0 Hs_dI mp

%% LLG parameters %%
%%%------------
q=1.6e-19;                              % Coulombs
hbar=6.626e-34/2/pi;                    % Reduced Planck's constant (J-s)
alpha = 0.1;                           %0.01 Gilbert damping parameter
g = 1.76e7*4*pi/1000;                   %1.76e7 Gyromagnetic ratio [(rad)/(Oe.s)]
mu0=1.2566e-6;                          %Magnetic permeability
uB=9.27400949e-24;
T=300;
kB=1.38e-23;

%% Magnet Parameters (taken from experiment)

Ms = 1.2e6; % A/m, Saturation Magnetization [emu/cm^3] emu/cm^3 = 1 KA/m
Ku2 = 66240; % Uni. anisotropy constant [J/m^3]
Hk = 2*Ku2/mu0/Ms; % Switching field [A/m]
Hk;
%%%%%%%%%%%%%     VCMA[MKS]    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%72.864%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=1;
V_vcma=1.048;                        %1.05984 HVCMA=HK Voltage across the MTJ,[V]
ksi=100e-15;                    %VCMA constant [J/V.m];37e-10[erg/V.cm]
t_fl=1e-9;                      %thickness of freelayer
Area=50*50*10^-18;
Vol=t_fl*Area;                  % Volume [m^3]
gamma=g/(1+alpha^2);            % reduced factor
d_MgO=1.6e-9;                   % thickness of MgO barrier[m]
% J_H_conv = hbar*P/(Hk*mu0*q*t_fl*Ms);
Hs_dI=hbar*P/(2*q*mu0*t_fl*Ms*Area);
H_VCMA_dv=2*ksi/(mu0*Ms*d_MgO*t_fl); %*mz,[Oe]
H_VCMA=V_vcma*H_VCMA_dv;
H_VCMA;
delta=(mu0*(Hk-H_VCMA)*Ms*Vol)/(2*kB*T);
delta;
% I_c=(alpha*g*q/uB/P)*Ms*Vol*(Hk-H_VCMA);    %from zhao Recent progresses in STT-MRAM 2016
I_c=2*q*alpha*2*delta*kB*T/hbar/P;        %from datta 2017
I=-0*10^-6;                            %current

Hext=0;                         % external field
mp=[0 0 1];                     % pinned layer
theta0=sqrt(kB*T/mu0/Ms/Hk/Vol);  % initial angle
% theta0=0.3;
A=2*delta*kB*T*5;
B=A/2.8;
step=1e-2;

%% Initial conditions of the simulation %%
mz = cos(theta0); % Magnet slightly off easy axis
m_initial = [0 sqrt(1-mz^2) mz];

%%   LLG     %%
time1 = 36000; % duration ns
tot1 = time1/step;
delta_t=step*10^-9; % time step
[theta,phi,Time,Energy]=deal(zeros(1,tot1));
[m,m22,m23,m24,m32,m33,m34,m42,m43,m44]=deal(zeros(tot1,3));
m_positive=0;
m_negtive=0;
for i=1:1:tot1
    Time(i)=i*step;
    if(i<36000/step)
        if(rem(i,10/step)==0)            %VCMA annealing
            V_vcma=V_vcma-0.000015;
            H_VCMA=2*ksi*V_vcma/(mu0*Ms*d_MgO*t_fl);
            delta=(mu0*(Hk-H_VCMA)*Ms*Vol)/(2*kB*T);
            delta
            A=2*delta*kB*T*3;
            B=A/2.8;
        end
    end
% 
%     if(i<30000/step)
%         if(rem(i,10/step)==0)              %temperature annealing
%             T=T-0.1;
%             T
%         end
%     else
%         if(rem(i,10/step)==0)              %temperature annealing
%             T=T-0.0001;
%             T
%         end
%     end
 
    %discription of thermal noise
    for j=1:1:9
        r=0+1.*randn(1,3); % Normal distribution, expectation = 0, standard deviation = 1
        Hth(j,:)=(sqrt(2*kB*T*alpha/(mu0*Ms*gamma*Vol*delta_t)))*r; % thermal noise field
    end

    Hs=Hs_dI*(2*q*alpha/hbar)*(A+1.183*B);
    if i==1
        [m22(i,:),m23(i,:),m24(i,:),m32(i,:),m33(i,:),m34(i,:),m42(i,:),m43(i,:),m44(i,:)]=deal(m_initial);
    else
        I22=-(2*q*alpha/hbar)*(A/2+A+1.183*B+A/2*(m23_read+m24_read+m32_read+m42_read)+B/4*m33_read+0.433*B*m43_read);
        I23=-(2*q*alpha/hbar)*(A/2+A+1.366*B+A/2*(m22_read+m24_read+m33_read+m43_read)+B/4*(m32_read+m34_read)+0.433*B*(m42_read+m44_read));
        I24=-(2*q*alpha/hbar)*(A/2+A+1.183*B+A/2*(m22_read+m23_read+m34_read+m44_read)+B/4*m33_read+0.433*B*m43_read);
        I32=-(2*q*alpha/hbar)*(A/2+A+B+A/2*(m33_read+m34_read+m22_read+m42_read)+B/4*m23_read+B/4*m43_read);
        I33=-(2*q*alpha/hbar)*(A/2+A+B+A/2*(m32_read+m34_read+m23_read+m43_read)+B/4*(m22_read+m24_read+m42_read+m44_read));
        I34=-(2*q*alpha/hbar)*(A/2+A+B+A/2*(m32_read+m33_read+m24_read+m44_read)+B/4*(m23_read+m43_read));
        I42=-(2*q*alpha/hbar)*(A/2+A+1.183*B+A/2*(m43_read+m44_read+m22_read+m32_read)+B/4*m33_read+0.433*B*m23_read);
        I43=-(2*q*alpha/hbar)*(A/2+A+1.366*B+A/2*(m42_read+m44_read+m23_read+m33_read)+0.433*B*(m22_read+m24_read)+B/4*(m32_read+m34_read));
        I44=-(2*q*alpha/hbar)*(A/2+A+1.183*B+A/2*(m43_read+m43_read+m24_read+m34_read)+0.433*B*m23_read+B/4*m33_read);

        m22(i,:)=sLLG_step(I22, m22(i-1,:), Hth(1,:), V_vcma, delta_t);
        m23(i,:)=sLLG_step(I23, m23(i-1,:), Hth(2,:), V_vcma, delta_t);
        m24(i,:)=sLLG_step(I24, m24(i-1,:), Hth(3,:), V_vcma, delta_t);
        m32(i,:)=sLLG_step(I32, m32(i-1,:), Hth(4,:), V_vcma, delta_t);
        m33(i,:)=sLLG_step(I33, m33(i-1,:), Hth(5,:), V_vcma, delta_t);
        m34(i,:)=sLLG_step(I34, m34(i-1,:), Hth(6,:), V_vcma, delta_t);
        m42(i,:)=sLLG_step(I42, m42(i-1,:), Hth(7,:), V_vcma, delta_t);
        m43(i,:)=sLLG_step(I43, m43(i-1,:), Hth(8,:), V_vcma, delta_t);
        m44(i,:)=sLLG_step(I44, m44(i-1,:), Hth(9,:), V_vcma, delta_t);

    end
    m22_read=sign(m22(i,3));    m23_read=sign(m23(i,3));    m24_read=sign(m24(i,3));    m32_read=sign(m32(i,3));
    m33_read=sign(m33(i,3));    m34_read=sign(m34(i,3));    m42_read=sign(m42(i,3));    m43_read=sign(m43(i,3));    m44_read=sign(m44(i,3));

    x22=0.5*(m22(i,3)+1);x23=0.5*(m23(i,3)+1);x24=0.5*(m24(i,3)+1);x32=0.5*(m32(i,3)+1);x33=0.5*(m33(i,3)+1);x34=0.5*(m34(i,3)+1);
    x42=0.5*(m42(i,3)+1);x43=0.5*(m43(i,3)+1);x44=0.5*(m44(i,3)+1);

    m_out=[m22_read, m23_read, m24_read;m32_read,m33_read,m34_read;m42_read,m43_read,m44_read];
    Energy(i)=2.8*((1-x22-x23-x24)^2+(1-x32-x33-x34)^2+(1-x42-x43-x44)^2+(1-x22-x32-x42)^2+(1-x23-x33-x43)^2+(1-x24-x34-x44)^2)+...
            1*(x22+x32+x42+x24+x34+x44+x22*x33+x23*x34+x32*x23+x33*x24+x32*x43+x33*x44+x42*x33+x43*x34+1.732*(x22*x43+x23*x44+x42*x23+x43*x24));
    %count
    % if(m(i,3)>0)
    %     m_positive=m_positive+1;
    % else
    %     m_negtive=m_negtive+1;
    % end
end

%% plot
m_out
% figure;
% hold on;
% plot(Time,m(:,1),'k-','linewidth',1.0); % m_x
% plot(Time,m(:,2),'b-','linewidth',1.0); % m_y
% hold on;
% plot(Time,m22(:,3),'-','linewidth',1.0); % m_z
% plot(Time,m23(:,3),'-','linewidth',1.0); % m_z
% plot(Time,m24(:,3),'-','linewidth',1.0); % m_z
% figure;
% plot(Time,m23(:,3),'-','linewidth',1.0); % m_z
% figure;
% plot(Time,m24(:,3),'-','linewidth',1.0); % m_z

figure;
plot(Time,Energy);title('energy1');
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%