%%%%   stochatic computing based neural network  benchmark iris   %%%
%%%%   writen by zuodong zhang  2018/3/26                          %%%
clear;
tic;
%% 参数定义
for x=1:1:50
    bits_length=64;         %随机数产生的精度
    sequence_length=32;    %随机数序列的长度
    
    iris_data = csvread('iris_test.csv',1,0);
    layer1_weights=[1.91197,	0.54162,	-3.87953,	2.06347;
                    5.04209,	1.75114,	-3.37571,	-6.10828;
                    -7.41239,	-3.74421,	8.69002,	7.8777;
                    -9.90473,	-4.98845,	11.0349,	6.84323];
    
    layer1_biases=[2.57034,     1.28579,    -2.67758,   -1.77147];
    
    logits_weights=[7.62891,	5.37066,	-9.66293;
                    4.45099,	-0.0872815,	-5.22292;
                    -6.80494,	-6.12244,	11.0609;
                    -9.24437,	5.95156,	6.54122];
    
    logits_biases=[2.39294,-0.601043,-0.988045];
    
    layer1_weights=abs(layer1_weights./12);
    logits_weights=abs(logits_weights./12);
    layer1_biases=abs(layer1_biases./12);
    logits_biases=abs(logits_biases./12);
    
    %%
    accuracy(x)=0;
    for i=1:1:30
        %     i=1;
        test=0.1*iris_data(i,1:4);
        [output_1,output_2,output_3]=deal(zeros(1,sequence_length));
        max_state=2*5-1;
        [state1,state2,state3,state4]=deal(floor(max_state/2));
        for j=1:1:sequence_length
            for m=1:1:4
                test_s(m)=s_bit(test(m)*bits_length,bits_length);
                layer1_biases_s(m)=s_bit(layer1_biases(m)*bits_length,bits_length);
            end
            for m=1:1:3
                logits_biases_s(m)=s_bit(logits_biases(m)*bits_length,bits_length);
            end
            
            
            for m=1:1:4
                for n=1:1:4
                    if(layer1_weights(m,n)<1)
                        layer1_weights_s(m,n)=s_bit(layer1_weights(m,n)*bits_length,bits_length);
                    else
                        layer1_weights_s(m,n)=s_bit((layer1_weights(m,n)-1)*bits_length,bits_length);
                    end
                end
            end
            for m=1:1:4
                for n=1:1:3
                    logits_weights_s(m,n)=s_bit(logits_weights(m,n)*bits_length,bits_length);
                end
            end
            %放缩15倍的计算
            voltage_1=(test_s(1)*layer1_weights_s(1,1)+test_s(2)*layer1_weights_s(2,1)-test_s(3)*layer1_weights_s(3,1)-test_s(4)*layer1_weights_s(4,1)+layer1_biases_s(1));
            voltage_2=(test_s(1)*layer1_weights_s(1,2)+test_s(2)*layer1_weights_s(2,2)-test_s(3)*layer1_weights_s(3,2)-test_s(4)*layer1_weights_s(4,2)+layer1_biases_s(2));
            voltage_3=(-test_s(1)*layer1_weights_s(1,3)-test_s(2)*layer1_weights_s(2,3)+test_s(3)*layer1_weights_s(3,3)+test_s(4)*layer1_weights_s(4,3)-layer1_biases_s(3));
            voltage_4=(test_s(1)*layer1_weights_s(1,4)-test_s(2)*layer1_weights_s(2,4)+test_s(3)*layer1_weights_s(3,4)+test_s(4)*layer1_weights_s(4,4)-layer1_biases_s(4));
            
            state1=state1+voltage_1;
            state2=state2+voltage_2;
            state3=state3+voltage_3;
            state4=state4+voltage_4;
            
            if(state1>max_state)                         state1=max_state;            end
            if(state1<0)                                    state1=0;            end
            if(state1>floor(max_state/2))              
                layer1_1=1;           
            else
                layer1_1=0;    
            end
            
            if(state2>max_state)                state2=max_state;            end
            if(state2<0)                state2=0;            end
            if(state2>floor(max_state/2))
                layer1_2=1;
            else
                layer1_2=0;
            end
            
            if(state3>max_state)                state3=max_state;            end
            if(state3<0)                state3=0;            end
            if(state3>floor(max_state/2))
                layer1_3=1;
            else
                layer1_3=0;
            end
            
            if(state4>max_state)                state4=max_state;            end
            if(state4<0)                state4=0;            end
            if(state4>floor(max_state/2))
                layer1_4=1;
            else
                layer1_4=0;
            end
            
            
            voltage2_1=(layer1_1*logits_weights_s(1,1)+layer1_2*logits_weights_s(2,1)-layer1_3*logits_weights_s(3,1)-layer1_4*logits_weights_s(4,1)+logits_biases_s(1));
            voltage2_2=(layer1_1*logits_weights_s(1,2)-layer1_2*logits_weights_s(2,2)-layer1_3*logits_weights_s(3,2)+layer1_4*logits_weights_s(4,2)-logits_biases_s(2));
            voltage2_3=(-layer1_1*logits_weights_s(1,3)-layer1_2*logits_weights_s(2,3)+layer1_3*logits_weights_s(3,3)+layer1_4*logits_weights_s(4,3)-logits_biases_s(3));
            
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