%%%%   stochatic computing based neural network  benchmark iris   %%%
%%%%   writen by zuodong zhang  2018/3/26                          %%%
clear;
tic;
%% 参数定义
for x=1:1:50
    bit=9;
    sequence_length=round(2^bit);    %随机数序列的长度
    
    iris_data = csvread('iris_test.csv',1,0);
    layer1_weights=[-1.0034,	2.11911,	-1.62791,	0.231651,	0.61993,	-1.83723,	-1.98954,	1.88547,	1.09951,	0.601555;
                    -2.31297,	2.67501,	-3.63939,	4.03773,	-1.5156,	-1.44924,	-1.25797,	3.04525,	3.58217,	-3.83612;
                    4.28507,	-5.71766,	4.98724,	-5.11524,	2.43317,	4.29179,	3.66089,	-4.86268,	-4.67096,	5.99592;
                    5.86411,	-7.64295,	5.87542,	-5.39865,	-1.51316,	6.6307,	5.89033,	-8.00004,	-5.9575,	4.30958];
    
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
        max_state=8-1;
        [state1,state2,state3,state4,state5,state6,state7,state8,state9,state10]=deal(floor(max_state/2));
        
        test_s1=LFSR(1,bit,test(1));test_s2=LFSR(2,bit,test(2));test_s3=LFSR(3,bit,test(3));test_s4=LFSR(4,bit,test(4));
        layer1_biases_s1=LFSR(5,bit,layer1_biases(1));layer1_biases_s2=LFSR(5,bit,layer1_biases(2));layer1_biases_s3=LFSR(5,bit,layer1_biases(3));layer1_biases_s4=LFSR(5,bit,layer1_biases(4));
        layer1_biases_s5=LFSR(5,bit,layer1_biases(5));layer1_biases_s6=LFSR(5,bit,layer1_biases(6));layer1_biases_s7=LFSR(5,bit,layer1_biases(7));layer1_biases_s8=LFSR(5,bit,layer1_biases(8));
        layer1_biases_s9=LFSR(5,bit,layer1_biases(9));layer1_biases_s10=LFSR(5,bit,layer1_biases(10));
        
        logits_biases_s1=LFSR(5,bit,logits_biases(1));logits_biases_s2=LFSR(5,bit,logits_biases(2));logits_biases_s3=LFSR(5,bit,logits_biases(3));
        
        layer1_weights_s11=LFSR(6,bit,layer1_weights(1,1));layer1_weights_s12=LFSR(6,bit,layer1_weights(1,2));layer1_weights_s13=LFSR(6,bit,layer1_weights(1,3));layer1_weights_s14=LFSR(6,bit,layer1_weights(1,4));
        layer1_weights_s15=LFSR(6,bit,layer1_weights(1,5));layer1_weights_s16=LFSR(6,bit,layer1_weights(1,6));layer1_weights_s17=LFSR(6,bit,layer1_weights(1,7));layer1_weights_s18=LFSR(6,bit,layer1_weights(1,8));
        layer1_weights_s19=LFSR(6,bit,layer1_weights(1,9));layer1_weights_s110=LFSR(6,bit,layer1_weights(1,10));
        
        layer1_weights_s21=LFSR(7,bit,layer1_weights(2,1));layer1_weights_s22=LFSR(7,bit,layer1_weights(2,2));layer1_weights_s23=LFSR(7,bit,layer1_weights(2,3));layer1_weights_s24=LFSR(7,bit,layer1_weights(2,4));
        layer1_weights_s25=LFSR(7,bit,layer1_weights(2,5));layer1_weights_s26=LFSR(7,bit,layer1_weights(2,6));layer1_weights_s27=LFSR(7,bit,layer1_weights(2,7));layer1_weights_s28=LFSR(7,bit,layer1_weights(2,8));
        layer1_weights_s29=LFSR(7,bit,layer1_weights(2,9));layer1_weights_s210=LFSR(7,bit,layer1_weights(2,10));
        
        layer1_weights_s31=LFSR(8,bit,layer1_weights(3,1));layer1_weights_s32=LFSR(8,bit,layer1_weights(3,2));layer1_weights_s33=LFSR(8,bit,layer1_weights(3,3));layer1_weights_s34=LFSR(8,bit,layer1_weights(3,4));
        layer1_weights_s35=LFSR(8,bit,layer1_weights(3,5));layer1_weights_s36=LFSR(8,bit,layer1_weights(3,6));layer1_weights_s37=LFSR(8,bit,layer1_weights(3,7));layer1_weights_s38=LFSR(8,bit,layer1_weights(3,8));
        layer1_weights_s39=LFSR(8,bit,layer1_weights(3,9));layer1_weights_s310=LFSR(8,bit,layer1_weights(3,10));
        
        layer1_weights_s41=LFSR(9,bit,layer1_weights(4,1));layer1_weights_s42=LFSR(9,bit,layer1_weights(4,2));layer1_weights_s43=LFSR(9,bit,layer1_weights(4,3));layer1_weights_s44=LFSR(9,bit,layer1_weights(4,4));
        layer1_weights_s45=LFSR(9,bit,layer1_weights(4,5));layer1_weights_s46=LFSR(9,bit,layer1_weights(4,6));layer1_weights_s47=LFSR(9,bit,layer1_weights(4,7));layer1_weights_s48=LFSR(9,bit,layer1_weights(4,8));
        layer1_weights_s49=LFSR(9,bit,layer1_weights(4,9));layer1_weights_s410=LFSR(9,bit,layer1_weights(4,10));
        
        logits_weights_s11=LFSR(10,bit,logits_weights(1,1));logits_weights_s12=LFSR(10,bit,logits_weights(1,2));logits_weights_s13=LFSR(10,bit,logits_weights(1,3));
        logits_weights_s21=LFSR(11,bit,logits_weights(2,1));logits_weights_s22=LFSR(11,bit,logits_weights(2,2));logits_weights_s23=LFSR(11,bit,logits_weights(2,3));
        logits_weights_s31=LFSR(12,bit,logits_weights(3,1));logits_weights_s32=LFSR(12,bit,logits_weights(3,2));logits_weights_s33=LFSR(12,bit,logits_weights(3,3));
        logits_weights_s41=LFSR(13,bit,logits_weights(4,1));logits_weights_s42=LFSR(13,bit,logits_weights(4,2));logits_weights_s43=LFSR(13,bit,logits_weights(4,3));
        logits_weights_s51=LFSR(14,bit,logits_weights(5,1));logits_weights_s52=LFSR(14,bit,logits_weights(5,2));logits_weights_s53=LFSR(14,bit,logits_weights(5,3));
        logits_weights_s61=LFSR(15,bit,logits_weights(6,1));logits_weights_s62=LFSR(15,bit,logits_weights(6,2));logits_weights_s63=LFSR(15,bit,logits_weights(6,3));
        logits_weights_s71=LFSR(16,bit,logits_weights(7,1));logits_weights_s72=LFSR(16,bit,logits_weights(7,2));logits_weights_s73=LFSR(16,bit,logits_weights(7,3));
        logits_weights_s81=LFSR(17,bit,logits_weights(8,1));logits_weights_s82=LFSR(17,bit,logits_weights(8,2));logits_weights_s83=LFSR(17,bit,logits_weights(8,3));
        logits_weights_s91=LFSR(18,bit,logits_weights(9,1));logits_weights_s92=LFSR(18,bit,logits_weights(9,2));logits_weights_s93=LFSR(18,bit,logits_weights(9,3));
        logits_weights_s101=LFSR(19,bit,logits_weights(10,1));logits_weights_s102=LFSR(19,bit,logits_weights(10,2));logits_weights_s103=LFSR(19,bit,logits_weights(10,3));
        
        for j=1:1:sequence_length
            
            voltage_1=(-test_s1(j)*layer1_weights_s11(j)-test_s2(j)*layer1_weights_s21(j)+test_s3(j)*layer1_weights_s31(j)+test_s4(j)*layer1_weights_s41(j)-layer1_biases_s1(j));
            voltage_2=(test_s1(j)*layer1_weights_s12(j)+test_s2(j)*layer1_weights_s22(j)-test_s3(j)*layer1_weights_s32(j)-test_s4(j)*layer1_weights_s42(j)+layer1_biases_s2(j));
            voltage_3=(-test_s1(j)*layer1_weights_s13(j)-test_s2(j)*layer1_weights_s23(j)+test_s3(j)*layer1_weights_s33(j)+test_s4(j)*layer1_weights_s43(j)-layer1_biases_s3(j));
            voltage_4=(test_s1(j)*layer1_weights_s14(j)+test_s2(j)*layer1_weights_s24(j)-test_s3(j)*layer1_weights_s34(j)-test_s4(j)*layer1_weights_s44(j)+layer1_biases_s4(j));
            voltage_5=(test_s1(j)*layer1_weights_s15(j)-test_s2(j)*layer1_weights_s25(j)+test_s3(j)*layer1_weights_s35(j)-test_s4(j)*layer1_weights_s45(j)-layer1_biases_s5(j));
            voltage_6=(-test_s1(j)*layer1_weights_s16(j)-test_s2(j)*layer1_weights_s26(j)+test_s3(j)*layer1_weights_s36(j)+test_s4(j)*layer1_weights_s46(j)-layer1_biases_s6(j));
            voltage_7=(-test_s1(j)*layer1_weights_s17(j)-test_s2(j)*layer1_weights_s27(j)+test_s3(j)*layer1_weights_s37(j)+test_s4(j)*layer1_weights_s47(j)-layer1_biases_s7(j));
            voltage_8=(test_s1(j)*layer1_weights_s18(j)+test_s2(j)*layer1_weights_s28(j)-test_s3(j)*layer1_weights_s38(j)-test_s4(j)*layer1_weights_s48(j)+layer1_biases_s8(j));
            voltage_9=(test_s1(j)*layer1_weights_s19(j)+test_s2(j)*layer1_weights_s29(j)-test_s3(j)*layer1_weights_s39(j)-test_s4(j)*layer1_weights_s49(j)+layer1_biases_s9(j));
            voltage_10=(test_s1(j)*layer1_weights_s110(j)-test_s2(j)*layer1_weights_s210(j)+test_s3(j)*layer1_weights_s310(j)+test_s4(j)*layer1_weights_s410(j)-layer1_biases_s10(j));

          %% FSM
            state1=state1+voltage_1;
            state2=state2+voltage_2;
            state3=state3+voltage_3;
            state4=state4+voltage_4;
            state5=state5+voltage_5;
            state6=state6+voltage_6;
            state7=state7+voltage_7;
            state8=state8+voltage_8;
            state9=state9+voltage_9;
            state10=state10+voltage_10;
            
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
            
            if(state5>max_state)                state5=max_state;            end
            if(state5<0)                state5=0;            end
            if(state5>floor(max_state/2))
                layer1_5=1;
            else
                layer1_5=0;
            end
            
            if(state6>max_state)                state6=max_state;            end
            if(state6<0)                state6=0;            end
            if(state6>floor(max_state/2))
                layer1_6=1;
            else
                layer1_6=0;
            end
            
            if(state7>max_state)                state7=max_state;            end
            if(state7<0)                state7=0;            end
            if(state7>floor(max_state/2))
                layer1_7=1;
            else
                layer1_7=0;
            end
            
            if(state8>max_state)                state8=max_state;            end
            if(state8<0)                state8=0;            end
            if(state8>floor(max_state/2))
                layer1_8=1;
            else
                layer1_8=0;
            end
            if(state9>max_state)                state9=max_state;            end
            if(state9<0)                state9=0;            end
            if(state9>floor(max_state/2))
                layer1_9=1;
            else
                layer1_9=0;
            end
            
            if(state10>max_state)                state10=max_state;            end
            if(state10<0)                state10=0;            end
            if(state10>floor(max_state/2))
                layer1_10=1;
            else
                layer1_10=0;
            end
            %%
            voltage2_1=(-layer1_1*logits_weights_s11(j)+layer1_2*logits_weights_s21(j)-layer1_3*logits_weights_s31(j)+layer1_4*logits_weights_s41(j)-layer1_5*logits_weights_s51(j)-layer1_6*logits_weights_s61(j)...
                -layer1_7*logits_weights_s71(j)+layer1_8*logits_weights_s81(j)+layer1_9*logits_weights_s91(j)-layer1_10*logits_weights_s101(j)-logits_biases_s1(j));
            voltage2_2=(-layer1_1*logits_weights_s12(j)+layer1_2*logits_weights_s22(j)-layer1_3*logits_weights_s32(j)-layer1_4*logits_weights_s42(j)+layer1_5*logits_weights_s52(j)-layer1_6*logits_weights_s62(j)...
                -layer1_7*logits_weights_s72(j)+layer1_8*logits_weights_s82(j)+layer1_9*logits_weights_s92(j)+layer1_10*logits_weights_s102(j)+logits_biases_s2(j));
            voltage2_3=(layer1_1*logits_weights_s13(j)-layer1_2*logits_weights_s23(j)+layer1_3*logits_weights_s33(j)-layer1_4*logits_weights_s43(j)-layer1_5*logits_weights_s53(j)+layer1_6*logits_weights_s63(j)...
                +layer1_7*logits_weights_s73(j)-layer1_8*logits_weights_s83(j)-layer1_9*logits_weights_s93(j)+layer1_10*logits_weights_s103(j)-logits_biases_s3(j));
            
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