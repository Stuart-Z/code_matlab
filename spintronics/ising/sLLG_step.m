function m= sLLG_step( I , m_pre ,Hth, V_vcma, delta_t)

global Hs_dI H_VCMA_dv Hk gamma alpha mp

Hs=Hs_dI*I;
Hthx=Hth(1);
Hthy=Hth(2);
Hthz=Hth(3);
H_VCMA=V_vcma*H_VCMA_dv;
H_eff=[ Hthx, Hthy, Hthz-H_VCMA*m_pre(3)+Hk*m_pre(3)];

m_next=m_pre+delta_t*gamma*(-cross(m_pre,H_eff)-alpha*cross(m_pre,cross(m_pre,H_eff))+...
    Hs*cross(cross(m_pre,mp),m_pre)+...
    alpha*Hs*cross(m_pre,mp));
theta =real(acos(m_next(3)));
phi =real(atan(m_next(2)/m_next(1)));
if(m_next(1)<=0)
    phi=phi+pi;
end
if(theta ==0)
    theta =0.01;
elseif(theta ==pi)
    theta =pi-0.01;
end
m(1)=sin(theta )*cos(phi);
m(2)=sin(theta )*sin(phi);
m(3)=cos(theta );

end

