%%% Writen by Tianqi Gao (2016/12/09)%%%
clear;
tic;
global gamma  alpha H_VCMA Hext Hk Ms...
    H_bias kB T mu0 Hth Hs_dI

%% LLG parameters %%
%%% Constants
%%%------------
q=1.6e-19;                              % Coulombs
hbar=6.626e-34/2/pi;                    % Reduced Planck's constant (J-s)
alpha = 0.1;                           %0.01 Gilbert damping parameter
g = 1.76e7*4*pi/1000;                   %1.76e7 Gyromagnetic ratio [(rad)/(Oe.s)]
mu0=1.2566e-6;                          %Magnetic permeability
T=300;
kB=1.38e-23;

%%% Magnet Parameters (taken from experiment)
%%%-------------------------------------------
Ms = 1.2e6; % A/m, Saturation Magnetization [emu/cm^3] emu/cm^3 = 1 KA/m
V = (50*50*1)*1e-27; % Volume [m^3]
Ku2 = 66240; % Uni. anisotropy constant [J/m^3]
Hk = 2*Ku2/mu0/Ms; % Switching field [A/m]

%%%%%%%%%%%%%     VCMA[MKS]    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%72.864%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=0.33;
Vs=1.0598;                        %0.4894��Voltage across the MTJ,[V]
ksi=100e-15;                    %VCMA constant [J/V.m];37e-10[erg/V.cm]
t_fl=1e-9;                      %thickness of freelayer
Area=50*50*10^-18;
Vol=t_fl*Area; % Volume [m^3]
gamma=g/(1+alpha^2);            % reduced factor
d_MgO=1.6e-9;                   % thickness of MgO barrier[m]
% J_H_conv = hbar*P/(Hk*mu0*q*t_fl*Ms);
Hs_dI=hbar*P/(2*q*mu0*t_fl*Ms*Area);
H_VCMA=2*ksi*Vs/(mu0*Ms*d_MgO*t_fl); %*mz,[Oe]
I=-0*10^-6;                            %current
Hext=0;                         % external field
H_bias=0;                       % bias magnetic field
mp=[0 0 1];                     % pinned layer
theta0=sqrt(kB*T/mu0/Ms/Hk/V);  % initial angle
% theta0=0.2;
step=1e-2;

%% Initial conditions of the simulation %%
mz = cos(theta0); % Magnet slightly off easy axis
m_initial = [0 sqrt(1-mz^2) mz];

%%   Step_1: Put the magnetization into plane   %%%%%%%
time1 = 10000; % step_1 duration ns
tot1 = time1/step;
delta_t=step*10^-9; % time step
[theta,phi,Time_pre]=deal(zeros(1,tot1));
m=zeros(tot1,3);
for i=1:1:tot1
    
    %discription of thermal noise
    r=0+1.*randn(1,3); % Normal distribution, expectation = 0, standard deviation = 1
    Hth=(sqrt(2*kB*T*alpha/(mu0*Ms*gamma*Vol*delta_t)))*r; % thermal noise field
    
    Hthx=Hth(1);
    Hthy=Hth(2);
    Hthz=Hth(3);
    
    Hs=Hs_dI*I;
    
    %     options = odeset('RelTol',1e-8,'AbsTol',1e-9);
    %     [t0,m1]=ode113('RNG_VCMA_Step_1',[0 0.001*1e-9*tau_c],m,options);
    if i==1
        m(i,:)=m_initial;
    else
        H_eff=[ Hthx, Hthy, Hthz-H_VCMA*m(i-1,3)+Hk*m(i-1,3)];
        m_next=m(i-1,:)+delta_t*gamma*(-cross(m(i-1,:),H_eff)-alpha*cross(m(i-1,:),cross(m(i-1,:),H_eff))+...
                    Hs*cross(cross(m(i-1,:),mp),m(i-1,:))+alpha*Hs*cross(m(i-1,:),mp));
        theta(i)=real(acos(m_next(3)));
        phi(i)=real(atan(m_next(2)/m_next(1)));
        if(theta(i)==0)
            theta(i)=0.01;
        elseif(theta(i)==pi)
            theta(i)=pi-0.01;
        end
        m(i,1)=sin(theta(i))*cos(phi(i));
        m(i,2)=sin(theta(i))*sin(phi(i));
        m(i,3)=cos(theta(i));
    end
    Time_pre(i)=i/1000;
end

%% Step_2: Generating random number %%%%%%%%%
% m_roll=m_pre(end,:);
% time2=0;% Step_2 duration
% tot2=time2/1e-3;
% m_eva=[];
% Time_roll=[];
% for i=1:1:tot2
%
%     r=0+1.*randn(1,3);
%     delta_t=1e-12;%time step
%     Vol=50*50*1*1e-27;%[m^3]
%     Hth=(sqrt(2*kB*T*alpha/(mu0*Ms*gamma*Vol*delta_t)))*r;
%
%     options = odeset('RelTol',1e-8,'AbsTol',1e-9);
%     [t,x]= ode113('RNG_Random_Step_2',[0 0.001*1e-9*tau_c], m_roll,options);
%
%     m_eva=[m_eva;x];
%     Time_roll=[Time_roll;(t+(i-1))];
%     m_roll = x(end,:);
% end

% T_eva=Time_pre(end,:)+Time_roll;
% T_tot=[Time_pre;T_eva];
% Mz=m_pre;
%
% tau=T_tot/1000;
% output=Mz(end,3);
%% plot

figure;
hold on;
% plot(Time_pre,m(:,1),'k-','linewidth',1.0); % m_x
% plot(Time_pre,m(:,2),'b-','linewidth',1.0); % m_y
plot(Time_pre,m(:,3),'r-','linewidth',1.0); % m_z
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%