function phi = f2( n )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
global alpha;
global hs;
if n==0
    phi=pi/2;
else
    phi=f2(n-1)-1-alpha*hs;
end
if phi>=2*pi
    phi=phi-pi*2;
end
end

