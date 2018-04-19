%%%%   stochatic computing based neural network  benchmark iris   %%%
%%%%   writen by zuodong zhang  2018/3/26                          %%%
clear;
tic;
%% 参数定义
for x=1:1:1
    bit=6;
    sequence_length=round(2^bit*0.9);    %随机数序列的长度
    
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
        max_state=8-1;
        [state1,state2,state3,state4]=deal(floor(max_state/2));
        
        test_s1=LFSR(1,bit,test(1));test_s2=LFSR(2,bit,test(2));test_s3=LFSR(3,bit,test(3));test_s4=LFSR(4,bit,test(4));
        layer1_biases_s1=LFSR(5,bit,layer1_biases(1));layer1_biases_s2=LFSR(6,bit,layer1_biases(2));layer1_biases_s3=LFSR(7,bit,layer1_biases(3));layer1_biases_s4=LFSR(8,bit,layer1_biases(4));
        logits_biases_s1=LFSR(10,bit,logits_biases(1));logits_biases_s2=LFSR(11,bit,logits_biases(2));logits_biases_s3=LFSR(12,bit,logits_biases(3));
        
        layer1_weights_s11=LFSR(21,bit,layer1_weights(1,1));layer1_weights_s12=LFSR(22,bit,layer1_weights(1,2));layer1_weights_s13=LFSR(23,bit,layer1_weights(1,3));layer1_weights_s14=LFSR(24,bit,layer1_weights(1,4));
        layer1_weights_s21=LFSR(31,bit,layer1_weights(2,1));layer1_weights_s22=LFSR(32,bit,layer1_weights(2,2));layer1_weights_s23=LFSR(33,bit,layer1_weights(2,3));layer1_weights_s24=LFSR(34,bit,layer1_weights(2,4));
        layer1_weights_s31=LFSR(41,bit,layer1_weights(3,1));layer1_weights_s32=LFSR(42,bit,layer1_weights(3,2));layer1_weights_s33=LFSR(43,bit,layer1_weights(3,3));layer1_weights_s34=LFSR(44,bit,layer1_weights(3,4));
        layer1_weights_s41=LFSR(51,bit,layer1_weights(4,1));layer1_weights_s42=LFSR(52,bit,layer1_weights(4,2));layer1_weights_s43=LFSR(53,bit,layer1_weights(4,3));layer1_weights_s44=LFSR(54,bit,layer1_weights(4,4));
        
        logits_weights_s11=LFSR(25,bit,logits_weights(1,1));logits_weights_s12=LFSR(26,bit,logits_weights(1,2));logits_weights_s13=LFSR(27,bit,logits_weights(1,3));
        logits_weights_s21=LFSR(35,bit,logits_weights(2,1));logits_weights_s22=LFSR(36,bit,logits_weights(2,2));logits_weights_s23=LFSR(37,bit,logits_weights(2,3));
        logits_weights_s31=LFSR(45,bit,logits_weights(3,1));logits_weights_s32=LFSR(46,bit,logits_weights(3,2));logits_weights_s33=LFSR(47,bit,logits_weights(3,3));
        logits_weights_s41=LFSR(55,bit,logits_weights(4,1));logits_weights_s42=LFSR(56,bit,logits_weights(4,2));logits_weights_s43=LFSR(57,bit,logits_weights(4,3));
        
        for j=1:1:sequence_length
            
            voltage_1=(test_s1(j)*layer1_weights_s11(j)+test_s2(j)*layer1_weights_s21(j)-test_s3(j)*layer1_weights_s31(j)-test_s4(j)*layer1_weights_s41(j)+layer1_biases_s1(j));
            voltage_2=(test_s1(j)*layer1_weights_s12(j)+test_s2(j)*layer1_weights_s22(j)-test_s3(j)*layer1_weights_s32(j)-test_s4(j)*layer1_weights_s42(j)+layer1_biases_s2(j));
            voltage_3=(-test_s1(j)*layer1_weights_s13(j)-test_s2(j)*layer1_weights_s23(j)+test_s3(j)*layer1_weights_s33(j)+test_s4(j)*layer1_weights_s43(j)-layer1_biases_s3(j));
            voltage_4=(test_s1(j)*layer1_weights_s14(j)-test_s2(j)*layer1_weights_s24(j)+test_s3(j)*layer1_weights_s34(j)+test_s4(j)*layer1_weights_s44(j)-layer1_biases_s4(j));
            
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
            
            
            voltage2_1=(layer1_1*logits_weights_s11(j)+layer1_2*logits_weights_s21(j)-layer1_3*logits_weights_s31(j)-layer1_4*logits_weights_s41(j)+logits_biases_s1(j));
            voltage2_2=(layer1_1*logits_weights_s12(j)-layer1_2*logits_weights_s22(j)-layer1_3*logits_weights_s32(j)+layer1_4*logits_weights_s42(j)-logits_biases_s2(j));
            voltage2_3=(-layer1_1*logits_weights_s13(j)-layer1_2*logits_weights_s23(j)+layer1_3*logits_weights_s33(j)+layer1_4*logits_weights_s43(j)-logits_biases_s3(j));
            
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