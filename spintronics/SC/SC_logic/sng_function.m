function [ out ] = sng_function( I ,V_apply)
%%%% Writen by Zuodong Zhang (2017/1/24) %%%%%

if(V_apply>1)
    V_apply=1;
end
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

Ms = 1.2e6; % A/m, Saturation Magnetization [emu/cm^3] emu/cm^3 = 1 KA/m  1.2e6
Ku2 = 66240; % Uni. anisotropy constant [J/m^3]
Hk = 2*Ku2/mu0/Ms; % Switching field [A/m]
%%%%%%%%%%%%%     VCMA[MKS]    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%72.864%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=0.33;
V_vcma=V_apply;                        %1.05984 HVCMA=HK Voltage across the MTJ,[V]
ksi=100e-15;                    %VCMA constant [J/V.m];37e-10[erg/V.cm]
t_fl=1e-9;                      %thickness of freelayer
Area=50*50*10^-18;
Vol=t_fl*Area; % Volume [m^3]
gamma=g/(1+alpha^2);            % reduced factor
d_MgO=1.6e-9;                   % thickness of MgO barrier[m]
% J_H_conv = hbar*P/(Hk*mu0*q*t_fl*Ms);
Hs_dI=hbar*P/(2*q*mu0*t_fl*Ms*Area);
H_VCMA_dv=2*ksi/(mu0*Ms*d_MgO*t_fl); %*mz,[Oe]
H_VCMA=V_vcma*H_VCMA_dv;
delta=(mu0*(Hk-H_VCMA)*Ms*Vol)/(2*kB*T);
% I_c=(alpha*g*q/uB/P)*Ms*Vol*(Hk-H_VCMA);    %from zhao Recent progresses in STT-MRAM 2016
I_c=2*q*alpha*2*delta*kB*T/hbar/P;        %from datta 2017
I=-I*10^-6;                            %current

Hext=0;                         % external field
mp=[0 0 1];                     % pinned layer
theta0=sqrt(kB*T/mu0/Ms/Hk/Vol);  % initial angle
% theta0=0.3;
step=1e-2;

%% Initial conditions of the simulation %%
mz = cos(theta0); % Magnet slightly off easy axis
m_initial = [0 sqrt(1-mz^2) mz];

%%   LLG     %%
time1 = 15;                 % step_1 duration ns
tot1 = time1/step;
delta_t=step*10^-9;         % time step
[theta,phi,Time]=deal(zeros(1,tot1));
m=zeros(tot1,3);
for i=1:1:tot1
    Time(i)=i*step;

    if(i>10/step)
        I=0;
    end
    %discription of thermal noise
    r=0+1.*randn(1,3);                                      % Normal distribution, expectation = 0, standard deviation = 1
    Hth=(sqrt(2*kB*T*alpha/(mu0*Ms*gamma*Vol*delta_t)))*r;  % thermal noise field
    
    Hthx=Hth(1);
    Hthy=Hth(2);
    Hthz=Hth(3);
    
    Hs=Hs_dI*I;
    
    if i==1
        m(i,:)=m_initial;
    else
        H_eff=[ Hthx, Hthy, Hthz-H_VCMA*m(i-1,3)+Hk*m(i-1,3)+Hext];
        m_next=m(i-1,:)+delta_t*gamma*(-cross(m(i-1,:),H_eff)-alpha*cross(m(i-1,:),cross(m(i-1,:),H_eff))+...
                    Hs*cross(cross(m(i-1,:),mp),m(i-1,:))+alpha*Hs*cross(m(i-1,:),mp));
        theta(i)=real(acos(m_next(3)));
        phi(i)=real(atan(m_next(2)/m_next(1)));
        if(m_next(1)<=0)
            phi(i)=phi(i)+pi;
        end
        if(theta(i)==0)
            theta(i)=0.01;
        elseif(theta(i)==pi)
            theta(i)=pi-0.01;
        end
        m(i,1)=sin(theta(i))*cos(phi(i));
        m(i,2)=sin(theta(i))*sin(phi(i));
        m(i,3)=cos(theta(i));
    end
end

%% plot
if(m(end,3)>0)
    out=0;
else
    out=1;
end

% figure;
% hold on;
% % plot(Time,m(:,1),'k-','linewidth',1.0); % m_x
% % plot(Time,m(:,2),'b-','linewidth',1.0); % m_y
% plot(Time,m(:,3),'r-','linewidth',1.0); % m_z

end

