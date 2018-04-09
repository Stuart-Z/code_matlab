function   theta = f1( n )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
global alpha;
global hs;
if n==1
    theta=0.01;
else
    theta=(1-alpha-hs)*f1(n-1);
end
if theta>pi
    theta=pi;
end

