clear;
tic;
bit_length=128;
sequence_length=1024;
for j=1:1:400
    %%
    input1=unifrnd(0,bit_length);    sympol1=sign(rand(1,1)-0.5);
    input2=unifrnd(0,bit_length);    sympol2=sign(rand(1,1)-0.5);
    input3=unifrnd(0,bit_length);    sympol3=sign(rand(1,1)-0.5);
    input4=unifrnd(0,bit_length);    sympol4=sign(rand(1,1)-0.5);
    input5=unifrnd(0,bit_length);    sympol5=sign(rand(1,1)-0.5);
    input6=unifrnd(0,bit_length);    sympol6=sign(rand(1,1)-0.5);
    %     input7=unifrnd(0,bit_length);    sympol7=sign(rand(1,1)-0.5);
    %     input8=unifrnd(0,bit_length);    sympol8=sign(rand(1,1)-0.5);
    %     input9=unifrnd(0,bit_length);    sympol9=sign(rand(1,1)-0.5);
    %     input10=unifrnd(0,bit_length);   sympol10=sign(rand(1,1)-0.5);
    %     input11=unifrnd(0,bit_length);   sympol11=sign(rand(1,1)-0.5);
    %     input12=unifrnd(0,bit_length);   sympol12=sign(rand(1,1)-0.5);
    %     input13=unifrnd(0,bit_length);   sympol13=sign(rand(1,1)-0.5);
    
    %% generate bit stream
    %     for i=1:1:sequence_length
    %         rnb1(i)=s_bit(input1,bit_length);
    %         rnb2(i)=s_bit(input2,bit_length);
    %         rnb3(i)=s_bit(input3,bit_length);
    %         rnb4(i)=s_bit(input4,bit_length);
    %         rnb5(i)=s_bit(input5,bit_length);
    % %         rnb6(i)=s_bit(input6,bit_length);
    % %         rnb7(i)=s_bit(input7,bit_length);
    % %         rnb8(i)=s_bit(input8,bit_length);
    % %         rnb9(i)=s_bit(input9,bit_length);
    % %         rnb10(i)=s_bit(input10,bit_length);
    % %         rnb11(i)=s_bit(input11,bit_length);
    % %         rnb12(i)=s_bit(input12,bit_length);
    % %         rnb13(i)=s_bit(input13,bit_length);
    %
    % %         voltage_c=(sympol1*rnb1(i)+sympol2*rnb2(i)+sympol3*rnb3(i)+sympol4*rnb4(i)+sympol5*rnb5(i)...
    % %             +sympol6*rnb6(i)+sympol7*rnb7(i)+sympol8*rnb8(i)+sympol9*rnb9(i)+sympol10*rnb10(i)+sympol11*rnb11(i)+sympol12*rnb12(i)+sympol13*rnb13(i))/22;
    % %         voltage_c=(sympol1*rnb1(i)+sympol2*rnb2(i)+sympol3*rnb3(i)+sympol4*rnb4(i)+sympol5*rnb5(i)+sympol6*rnb6(i)+sympol7*rnb7(i)+sympol8*rnb8(i)+sympol9*rnb9(i))/6;
    % %         voltage_c=(sympol1*rnb1(i)+sympol2*rnb2(i)+sympol3*rnb3(i)+sympol4*rnb4(i)+sympol5*rnb5(i)+sympol6*rnb6(i)+sympol7*rnb7(i))/4;
    %         voltage_c=(sympol1*rnb1(i)+sympol2*rnb2(i)+sympol3*rnb3(i)+sympol4*rnb4(i)+sympol5*rnb5(i))/5;
    %
    % %         voltage_c=(sympol1*rnb1(i)+sympol2*rnb2(i))/3;
    % %         out(i)=sng_function_nn(340,voltage_c);
    %         out(i)=sigmf(voltage_c,[6,0]);
    %     end
    
    max_state=2*5-1;
    state=floor(max_state/2);
    %% FSM
    for i=1:1:sequence_length
        rnb1(i)=s_bit(input1,bit_length);
        rnb2(i)=s_bit(input2,bit_length);
        rnb3(i)=s_bit(input3,bit_length);
        rnb4(i)=s_bit(input4,bit_length);
        rnb5(i)=s_bit(input5,bit_length);
        rnb6(i)=s_bit(input6,bit_length);
        %         rnb7(i)=s_bit(input7,bit_length);
        %         rnb8(i)=s_bit(input8,bit_length);
        %         rnb9(i)=s_bit(input9,bit_length);
        %         rnb10(i)=s_bit(input10,bit_length);
        %         rnb11(i)=s_bit(input11,bit_length);
        %         rnb12(i)=s_bit(input12,bit_length);
        %         rnb13(i)=s_bit(input13,bit_length);
        
        sum_number=(sympol1*rnb1(i)+sympol2*rnb2(i)+sympol3*rnb3(i)+sympol4*rnb4(i)+sympol5*rnb5(i));
        state=state+sum_number;
        if(state>max_state)
            state=max_state;
        end
        if(state<0)
            state=0;
        end
        if(state>floor(max_state/2))
            out(i)=1;
        else
            out(i)=0;
        end
    end
    
    
    %% count
    output=sum(out);
    %     input(j)=(sympol1*input1+sympol2*input2+sympol3*input3+sympol4*input4+sympol5*input5...
    %         +sympol6*input6+sympol7*input7+sympol8*input8+sympol9*input9+sympol10*input10+sympol11*input11+sympol12*input12+sympol13*input13)/bit_length;
    %     input(j)=(sympol1*input1+sympol2*input2+sympol3*input3+sympol4*input4+sympol5*input5+sympol6*input6+sympol7*input7+sympol8*input8+sympol9*input9)/bit_length;
    %     input(j)=(sympol1*input1+sympol2*input2+sympol3*input3+sympol4*input4+sympol5*input5+sympol6*input6+sympol7*input7)/bit_length;
%     input(j)=(sympol1*input1+sympol2*input2+sympol3*input3+sympol4*input4+sympol5*input5+sympol6*input6)/bit_length;
    %     input(j)=(sympol1*input1+sympol2*input2)/bit_length;
    input(j)=(sympol1*input1+sympol2*input2+sympol3*input3+sympol4*input4+sympol5*input5)/bit_length;
    result(j)=output/sequence_length;
    sigmoid(j)=sigmf(input(j),[12,0]);
%     tan(j)=tanh(2*input(j));
    j
end
figure;
plot(input,result,'ro','MarkerSize',4);

hold on;
plot(input,sigmoid,'go','MarkerSize',4);
toc