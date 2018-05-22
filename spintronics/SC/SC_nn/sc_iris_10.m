%%%%   stochatic computing based neural network  benchmark iris   %%%
%%%%   writen by zuodong zhang  2018/3/26                          %%%
clear;
tic;
%% 参数定义
for x=1:1:50
    bits_length=64;         %随机数产生的精度
    sequence_length=64;    %随机数序列的长度
    
    iris_data = csvread('iris_test.csv',1,0);
    layer1_weights=[-1.0034,	2.11911,	-1.62791,	0.231651,	0.61993,	-1.83723,	-1.98954,	1.88547,	1.09951,	0.601555;
                    -2.31297,	2.67501,	-3.63939,	4.03773,	-1.5156,	-1.44924,	-1.25797,	3.04525,	3.58217,	-3.83612;
                    4.28507,	-5.71766,	4.98724,	-5.11524,	2.43317,	4.29179,	3.66089,	-4.86268,	-4.67096,	5.99592;
                    5.86411,	-7.64295,	5.87542,	-5.39865,	-1.51316,	6.6307,     5.89033,	-8.00004,	-5.9575,	4.30958];
    
    layer1_biases=[-1.55542,    1.99542,    -0.901998,  1.20303,    -0.485448,  -1.6638,    -1.14865,   1.66095,    1.18035,    -1.14597];
    
    logits_weights=[-4.95559,	-1.64442,	5.20687;
                    5.25891,	4.47529,	-7.07294;
                    -5.29129,	-0.172276,	6.47114;
                    6.06033,	-0.595987,	-5.31683;
                    -1.50574,	2.67004,	-0.467336;
                    -3.62322,	-3.19926,	5.57017;
                    -3.18659,	-2.16822,	5.83613;
                    5.02217,	3.75108,	-7.19316;
                    5.50535,	0.708707,	-6.07905;
                    -5.54193,	3.4235,     4.801];
    
    logits_biases=[-0.113811,   0.397479,   -0.1783];
    
    layer1_weights=abs(layer1_weights./8);
    logits_weights=abs(logits_weights./8);
    layer1_biases=abs(layer1_biases./8);
    logits_biases=abs(logits_biases./8);
    
    %%
    accuracy(x)=0;
    for i=1:1:30
        test=0.1*iris_data(i,1:4);
        [voltage_1,voltage_2,voltage_3,voltage_4,voltage_5,voltage_6,voltage_7,voltage_8,voltage_9,voltage_10]=deal(0);
        delta_v=0.2;
        for j=1:1:sequence_length
            for m=1:1:4
                                test_s(m)=s_bit(test(m)*bits_length,bits_length);
                                layer1_biases_s(m)=s_bit(layer1_biases(m)*bits_length,bits_length);
%                 test_s(m)=s_bit_MTJ(test(m));
%                 layer1_biases_s(m)=s_bit_MTJ(layer1_biases(m));
            end
            for m=1:1:10
                  layer1_biases_s(m)=s_bit(layer1_biases(m)*bits_length,bits_length);
%                 layer1_biases_s(m)=s_bit_MTJ(layer1_biases(m));
            end
            for m=1:1:3
                                logits_biases_s(m)=s_bit(logits_biases(m)*bits_length,bits_length);
%                 logits_biases_s(m)=s_bit_MTJ(logits_biases(m));
            end
            
            
            for m=1:1:4
                for n=1:1:10
                    if(layer1_weights(m,n)<1)
                                                layer1_weights_s(m,n)=s_bit(layer1_weights(m,n)*bits_length,bits_length);
%                         layer1_weights_s(m,n)=s_bit_MTJ(layer1_weights(m,n));
                    else
                                                layer1_weights_s(m,n)=s_bit((layer1_weights(m,n)-1)*bits_length,bits_length);
%                         layer1_weights_s(m,n)=s_bit_MTJ(layer1_weights(m,n)-1);
                    end
                end
            end
            for m=1:1:10
                for n=1:1:3
                                        logits_weights_s(m,n)=s_bit(logits_weights(m,n)*bits_length,bits_length);
%                     logits_weights_s(m,n)=s_bit_MTJ(logits_weights(m,n));
                end
            end
            
            voltage_1=(-test_s(1)*layer1_weights_s(1,1)-test_s(2)*layer1_weights_s(2,1)+test_s(3)*layer1_weights_s(3,1)+test_s(4)*layer1_weights_s(4,1)-layer1_biases_s(1))*delta_v+voltage_1;
            voltage_2=(test_s(1)*layer1_weights_s(1,2)+test_s(2)*layer1_weights_s(2,2)-test_s(3)*layer1_weights_s(3,2)-test_s(4)*layer1_weights_s(4,2)+layer1_biases_s(2))*delta_v+voltage_2;
            voltage_3=(-test_s(1)*layer1_weights_s(1,3)-test_s(2)*layer1_weights_s(2,3)+test_s(3)*layer1_weights_s(3,3)+test_s(4)*layer1_weights_s(4,3)-layer1_biases_s(3))*delta_v+voltage_3;
            voltage_4=(test_s(1)*layer1_weights_s(1,4)+test_s(2)*layer1_weights_s(2,4)-test_s(3)*layer1_weights_s(3,4)-test_s(4)*layer1_weights_s(4,4)+layer1_biases_s(4))*delta_v+voltage_4;
            voltage_5=(test_s(1)*layer1_weights_s(1,5)-test_s(2)*layer1_weights_s(2,5)+test_s(3)*layer1_weights_s(3,5)-test_s(4)*layer1_weights_s(4,5)-layer1_biases_s(5))*delta_v+voltage_5;
            voltage_6=(-test_s(1)*layer1_weights_s(1,6)-test_s(2)*layer1_weights_s(2,6)+test_s(3)*layer1_weights_s(3,6)+test_s(4)*layer1_weights_s(4,6)-layer1_biases_s(6))*delta_v+voltage_6;
            voltage_7=(-test_s(1)*layer1_weights_s(1,7)-test_s(2)*layer1_weights_s(2,7)+test_s(3)*layer1_weights_s(3,7)+test_s(4)*layer1_weights_s(4,7)-layer1_biases_s(7))*delta_v+voltage_7;
            voltage_8=(test_s(1)*layer1_weights_s(1,8)+test_s(2)*layer1_weights_s(2,8)-test_s(3)*layer1_weights_s(3,8)-test_s(4)*layer1_weights_s(4,8)+layer1_biases_s(8))*delta_v+voltage_8;
            voltage_9=(test_s(1)*layer1_weights_s(1,9)+test_s(2)*layer1_weights_s(2,9)-test_s(3)*layer1_weights_s(3,9)-test_s(4)*layer1_weights_s(4,9)+layer1_biases_s(9))*delta_v+voltage_9;
            voltage_10=(test_s(1)*layer1_weights_s(1,10)-test_s(2)*layer1_weights_s(2,10)+test_s(3)*layer1_weights_s(3,10)+test_s(4)*layer1_weights_s(4,10)-layer1_biases_s(10))*delta_v+voltage_10;
            
            if(voltage_1<-0.9) voltage_1=-0.9;  end
            if(voltage_1>0.9) voltage_1=0.9;    end
            if(voltage_2<-0.9) voltage_2=-0.9;  end
            if(voltage_2>0.9) voltage_2=0.9;    end
            if(voltage_3<-0.9) voltage_3=-0.9;  end
            if(voltage_3>0.9) voltage_3=0.9;    end
            if(voltage_4<-0.9) voltage_4=-0.9;  end
            if(voltage_4>0.9) voltage_4=0.9;    end
            if(voltage_5<-0.9) voltage_5=-0.9;  end
            if(voltage_5>0.9) voltage_5=0.9;    end
            if(voltage_6<-0.9) voltage_6=-0.9;  end
            if(voltage_6>0.9) voltage_6=0.9;    end
            if(voltage_7<-0.9) voltage_7=-0.9;  end
            if(voltage_7>0.9) voltage_7=0.9;    end
            if(voltage_8<-0.9) voltage_8=-0.9;  end
            if(voltage_8>0.9) voltage_8=0.9;    end
            if(voltage_9<-0.9) voltage_9=-0.9;  end
            if(voltage_9>0.9) voltage_9=0.9;    end
            if(voltage_10<-0.9) voltage_10=-0.9;  end
            if(voltage_10>0.9) voltage_10=0.9;    end

   
            
            layer1_1=s_bit(bits_length*sigmf(voltage_1,[6,0]),bits_length);
            layer1_2=s_bit(bits_length*sigmf(voltage_2,[6,0]),bits_length);
            layer1_3=s_bit(bits_length*sigmf(voltage_3,[6,0]),bits_length);
            layer1_4=s_bit(bits_length*sigmf(voltage_4,[6,0]),bits_length);
            layer1_5=s_bit(bits_length*sigmf(voltage_5,[6,0]),bits_length);
            layer1_6=s_bit(bits_length*sigmf(voltage_6,[6,0]),bits_length);
            layer1_7=s_bit(bits_length*sigmf(voltage_7,[6,0]),bits_length);
            layer1_8=s_bit(bits_length*sigmf(voltage_8,[6,0]),bits_length);
            layer1_9=s_bit(bits_length*sigmf(voltage_9,[6,0]),bits_length);
            layer1_10=s_bit(bits_length*sigmf(voltage_10,[6,0]),bits_length);
            
            %             layer1_1=sng_function_nn(340,voltage_1);
            %             layer1_2=sng_function_nn(340,voltage_2);
            %             layer1_3=sng_function_nn(340,voltage_3);
            %             layer1_4=sng_function_nn(340,voltage_4);
            
            voltage2_1=(-layer1_1*logits_weights_s(1,1)+layer1_2*logits_weights_s(2,1)-layer1_3*logits_weights_s(3,1)+layer1_4*logits_weights_s(4,1)-layer1_5*logits_weights_s(5,1)-layer1_6*logits_weights_s(6,1)...
                -layer1_7*logits_weights_s(7,1)+layer1_8*logits_weights_s(8,1)+layer1_9*logits_weights_s(9,1)-layer1_10*logits_weights_s(10,1)-logits_biases_s(1));
            voltage2_2=(-layer1_1*logits_weights_s(1,2)+layer1_2*logits_weights_s(2,2)-layer1_3*logits_weights_s(3,2)-layer1_4*logits_weights_s(4,2)+layer1_5*logits_weights_s(5,2)-layer1_6*logits_weights_s(6,2)...
                -layer1_7*logits_weights_s(7,2)+layer1_8*logits_weights_s(8,2)+layer1_9*logits_weights_s(9,2)+layer1_10*logits_weights_s(10,2)+logits_biases_s(2));
            voltage2_3=(layer1_1*logits_weights_s(1,3)-layer1_2*logits_weights_s(2,3)+layer1_3*logits_weights_s(3,3)-layer1_4*logits_weights_s(4,3)-layer1_5*logits_weights_s(5,3)+layer1_6*logits_weights_s(6,3)...
                +layer1_7*logits_weights_s(7,3)-layer1_8*logits_weights_s(8,3)-layer1_9*logits_weights_s(9,3)+layer1_10*logits_weights_s(10,3)-logits_biases_s(3));
            
            output_1(j)=voltage2_1;
            output_2(j)=voltage2_2;
            output_3(j)=voltage2_3;
            
        end
        
        output1=sum(output_1)/sequence_length;
        output2=sum(output_2)/sequence_length;
        output3=sum(output_3)/sequence_length;
        output=[output1,output2,output3];
        
        
        [temp,index]=max(output);
        if (index==iris_data(i,5)+1)
            accuracy(x)=accuracy(x)+1;
        end
        
    end
    x
end
mean(accuracy)
figure;
toc