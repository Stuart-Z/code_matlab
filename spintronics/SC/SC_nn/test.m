clear;
tic;
input=3;
bit_length=256;
sequence=1024;
for i=1:1:sequence
    in(i)=s_bit_MTJ(input/bit_length);
    
end
result=sum(in)/sequence*bit_length;
toc