function sequence = LFSR( init, bit_length,input)

n = bit_length;
N = 2^n;
register = dec2bin(init,bit_length)-'0';                     %定义移位寄存器的初始状态
register=zeros(1,n)+register;
newregister = zeros(1, n);
sequence = zeros(1, N);
switch(bit_length)
    case 10
        feedback=[0,0,1,0,0,0,0,0,0,1];
    case 9
        feedback=[0,0,0,1,0,0,0,0,1]; 
    case 8
        feedback=[0,1,1,1,0,0,0,1]; 
    case 7
        feedback=[0,0,1,0,0,0,1];  
    case 6
        feedback=[1,0,0,0,0,1];
    case 5
        feedback=[0,1,0,0,1];
end


for i = 1:N
%     compare(i,:)=register;
    register_num=0;
    for j=1:1:bit_length
        register_num=register_num+register(j)*2^(bit_length-j);
    end
    if(input*N>register_num)
        sequence(i)=1;
    else
        sequence(i)=0;
    end
    newregister(1)= mod(sum(feedback.*register),2);
    newregister(2:n) = register(1:n-1);
    register = newregister;
end

end

