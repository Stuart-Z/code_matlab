clear;
%% TSP problem simplify
% syms x11 x12 x13 x14 x21 x22 x23 x24 x31 x32 x33 x34 x41 x42 x43 x44 A B;
% syms m11 m12 m13 m14 m21 m22 m23 m24 m31 m32 m33 m34 m41 m42 m43 m44
% h=A*((1-x22-x32-x42)^2+(1-x23-x33-x43)^2+(1-x24-x34-x44)^2+(1-x22-x23-x24)^2+(1-x32-x33-x34)^2+(1-x42-x43-x44)^2)+...
%     B*(x22+x32+x42+x24+x34+x44+x22*x33+x23*x34+x32*x23+x33*x24+x32*x43+x33*x44+x42*x33+x43*x34+1.732*(x22*x43+x23*x44+x42*x23+x43*x24));
% h1=subs(h,{x22, x23, x24, x32, x33, x34, x42, x43, x44},{(m22+1)/2,(m23+1)/2,(m24+1)/2,(m32+1)/2,(m33+1)/2,(m34+1)/2,(m42+1)/2,(m43+1)/2,(m44+1)/2});
% h2=expand(h1)

%% energy 
for i=1:1:511
    x22=bitget(i,1);
    x23=bitget(i,2);
    x24=bitget(i,3);
    x32=bitget(i,4);
    x33=bitget(i,5);
    x34=bitget(i,6);
    x42=bitget(i,7);
    x43=bitget(i,8);
    x44=bitget(i,9);
    
    Energy(i)=2.8*((1-x22-x23-x24)^2+(1-x32-x33-x34)^2+(1-x42-x43-x44)^2+(1-x22-x32-x42)^2+(1-x23-x33-x43)^2+(1-x24-x34-x44)^2)+...
            1*(x22+x32+x42+x24+x34+x44+x22*x33+x23*x34+x32*x23+x33*x24+x32*x43+x33*x44+x42*x33+x43*x34+1.732*(x22*x43+x23*x44+x42*x23+x43*x24));
end
plot(Energy)