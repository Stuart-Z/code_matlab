%%%%   stochatic computing based neural network  benchmark iris   %%%
%%%%   writen by zuodong zhang  2018/3/26                          %%%
clear;
tic;
%% 参数定义
for x=1:1:1
bits_length=64;         %随机数产生的精度
sequence_length=128;    %随机数序列的长度

iris_data = csvread('iris_test.csv',1,0);
layer1_weights=[6.9032,	    5.04648,	-4.34227,	-5.02293;
                -6.24207,	-5.23362,	5.34296,	3.99471;
                14.8788,	13.9252,	-13.9821,	-14.523;
                11.3529,	12.4055,	-11.7516,	-12.2361];

layer1_biases=[-0.889334,   -0.844831,  0.799981,   0.926365];

logits_weights=[-13.7329,	1.27417,	14.1756;
                -13.4678,	-1.06278,	13.8347;
                12.981,     1.1352,     -13.7875;
                12.9109,	1.89396,	-13.8149];

logits_biases=[0.0456089,   0.852301,   -0.779679];

layer1_weights=abs(layer1_weights./15);
logits_weights=abs(logits_weights./15);
layer1_biases=abs(layer1_biases./1.5);
logits_biases=abs(logits_biases./15);

%%
accuracy(x)=0;
for i=1:1:30
%     i=1;
    test=0.1*iris_data(i,1:4);
    [output_1,output_2,output_3]=deal(zeros(1,sequence_length));
    for j=1:1:sequence_length
        for m=1:1:4
            test_s(m)=s_bit(test(m)*bits_length,bits_length);
            layer1_biases_s(m)=s_bit(layer1_biases(m)*bits_length,bits_length);
%             test_s(m)=s_bit_MTJ(test(m));
%             layer1_biases_s(m)=s_bit_MTJ(layer1_biases(m));
        end 
        for m=1:1:3
            logits_biases_s(m)=s_bit(logits_biases(m)*bits_length,bits_length);
%             logits_biases_s(m)=s_bit_MTJ(logits_biases(m));
        end


        for m=1:1:4
            for n=1:1:4
                if(layer1_weights(m,n)<1)
                    layer1_weights_s(m,n)=s_bit(layer1_weights(m,n)*bits_length,bits_length);
%                     layer1_weights_s(m,n)=s_bit_MTJ(layer1_weights(m,n));
                else
                    layer1_weights_s(m,n)=s_bit((layer1_weights(m,n)-1)*bits_length,bits_length);
%                     layer1_weights_s(m,n)=s_bit_MTJ(layer1_weights(m,n)-1);
                end
            end
        end
        for m=1:1:4
            for n=1:1:3
                logits_weights_s(m,n)=s_bit(logits_weights(m,n)*bits_length,bits_length);
%                 logits_weights_s(m,n)=s_bit_MTJ(logits_weights(m,n));
            end
        end

        voltage_1=(test_s(1)*layer1_weights_s(1,1)-test_s(2)*layer1_weights_s(2,1)+test_s(3)*layer1_weights_s(3,1)+test_s(4)*layer1_weights_s(4,1)-layer1_biases_s(1))/2;
        voltage_2=(test_s(1)*layer1_weights_s(1,2)-test_s(2)*layer1_weights_s(2,2)+test_s(3)*layer1_weights_s(3,2)+test_s(4)*layer1_weights_s(4,2)-layer1_biases_s(2))/2;
        voltage_3=-(test_s(1)*layer1_weights_s(1,3)-test_s(2)*layer1_weights_s(2,3)+test_s(3)*layer1_weights_s(3,3)+test_s(4)*layer1_weights_s(4,3)-layer1_biases_s(3))/2;
        voltage_4=-(test_s(1)*layer1_weights_s(1,4)-test_s(2)*layer1_weights_s(2,4)+test_s(3)*layer1_weights_s(3,4)+test_s(4)*layer1_weights_s(4,4)-layer1_biases_s(4))/2;

        
            layer1_1=s_bit(bits_length*sigmf(voltage_1,[6,0]),bits_length);
            layer1_2=s_bit(bits_length*sigmf(voltage_2,[6,0]),bits_length);
            layer1_3=s_bit(bits_length*sigmf(voltage_3,[6,0]),bits_length);
            layer1_4=s_bit(bits_length*sigmf(voltage_4,[6,0]),bits_length);
        
%         layer1_1=sng_function_nn(340,voltage_1);
%         layer1_2=sng_function_nn(340,voltage_2);
%         layer1_3=sng_function_nn(340,voltage_3);
%         layer1_4=sng_function_nn(340,voltage_4);

        voltage2_1=(-layer1_1*logits_weights_s(1,1)-layer1_2*logits_weights_s(2,1)+layer1_3*logits_weights_s(3,1)+layer1_4*logits_weights_s(4,1)+logits_biases_s(1))/2;
        voltage2_2=(layer1_1*logits_weights_s(1,2)-layer1_2*logits_weights_s(2,2)+layer1_3*logits_weights_s(3,2)+layer1_4*logits_weights_s(4,2)+logits_biases_s(2))/2;
        voltage2_3=(layer1_1*logits_weights_s(1,3)+layer1_2*logits_weights_s(2,3)-layer1_3*logits_weights_s(3,3)-layer1_4*logits_weights_s(4,3)-logits_biases_s(3))/2;

        output_1(j)=voltage2_1;
        output_2(j)=voltage2_2;
        output_3(j)=voltage2_3;

        % output_1=layer1_1*logits_weights(1,1)+layer1_2*logits_weights(2,1)+layer1_3*logits_weights(3,1)+layer1_4*logits_weights(4,1)+logits_biases(1);
        % output_2=layer1_1*logits_weights(1,2)+layer1_2*logits_weights(2,2)+layer1_3*logits_weights(3,2)+layer1_4*logits_weights(4,2)+logits_biases(2);
        % output_3=layer1_1*logits_weights(1,3)+layer1_2*logits_weights(2,3)+layer1_3*logits_weights(3,3)+layer1_4*logits_weights(4,3)+logits_biases(3);    
    end

    output1=sum(output_1)/sequence_length;
    output2=sum(output_2)/sequence_length;
    output3=sum(output_3)/sequence_length;
    output=[output1,output2,output3];
    
    
    [temp,index]=max(output);
    if (index==iris_data(i,5)+1)
        accuracy(x)=accuracy(x)+1;
    end
    i
end
x
end
mean(accuracy)
figure;
toc