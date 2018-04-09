function dmdt = RNG_Random_Step_2(t,m)
global J_H_conv alpha H_VCMA Hk Ms...
    H_bias Hth

%%% demagnetizing factor
Nx=0*0.0045;
Ny=0*0.0152;
Nz=0*1;%0.9803;
H_demx=-Ms*Nx*m(1)*4*pi/1000;
H_demy=-Ms*Ny*m(2)*4*pi/1000;
H_demz=-Ms*Nz*m(3)*4*pi/1000;

%%% Thermal field
Hthx=Hth(1);
Hthy=Hth(2);
Hthz=Hth(3);

%%% density of STT current
J=0;%A/m^2

%%% effective field
H_eff=[0*H_bias+H_demx+Hthx,H_demy+Hthy,Hthz+H_demz+...
    H_VCMA*m(3)+Hk*m(3)];
h_eff = H_eff/Hk;
%%% direction of fixed layer
mp = [0 0 1];

%%% Differential Equation for magnetization Dynamics
dmdt0=(1/(2*alpha))*(-cross(m,h_eff)-alpha*cross(m,cross(m,h_eff))+J_H_conv*J*cross(m,cross(m,mp)));
dmdt=dmdt0';
end