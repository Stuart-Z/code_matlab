function phi = f2( n )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
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

