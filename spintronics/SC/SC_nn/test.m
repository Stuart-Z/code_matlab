clear;
% tic;
% input=3;
% bit_length=256;
% sequence=1024;
% for i=1:1:sequence
%     in(i)=s_bit_MTJ(input/bit_length);
%
% end
% result=sum(in)/sequence*bit_length;
% toc

input=0.5;
a=LFSR(50,10,input);
sum(a)/length(a)
% bit_length=6;
% for k=1:1:2^6-1
%     init=k;
%     n = bit_length;
%     N = 2^n;
%     register = dec2bin(init,bit_length)-'0';                     %������λ�Ĵ����ĳ�ʼ״̬
%     register=zeros(1,n)+register;
%     newregister = zeros(1, n);
%     sequence = zeros(1, N);
% %     feedback=[0,0,1,0,0,0,0,0,0,1];                              %10bit�ķ���
% %     feedback=[0,0,0,1,0,0,0,0,1];                               %9bit�ķ���
% %     feedback=[0,1,1,1,0,0,0,1];                                 %8bit�ķ���
% %     feedback=[0,0,1,0,0,0,1];                                   %7bit�ķ���
%     feedback=[1,0,0,0,0,1];                                     %6bit�ķ���
% 
%     
%     for i = 1:N
%         compare(i,:)=register;
%         register_num=0;
%         for j=1:1:bit_length
%             register_num=register_num+register(j)*2^(bit_length-j);
%         end
%         if(input*N>register_num)
%             sequence(i)=1;
%         else
%             sequence(i)=0;
%         end
%         newregister(1)= mod(sum(feedback.*register),2);
%         newregister(2:n) = register(1:n-1);
%         register = newregister;
%     end
%     error(k)=sum(sequence)/N-input;
%     
% end
% plot(error);