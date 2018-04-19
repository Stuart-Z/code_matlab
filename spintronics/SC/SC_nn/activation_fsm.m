clear;
tic;
bit=10;
bit_length=2^bit;
sequence_length=1024;
for j=1:1:200
    %%
    input1=unifrnd(0,bit_length)/bit_length;    sympol1=sign(rand(1,1)-0.5);
    input2=unifrnd(0,bit_length)/bit_length;    sympol2=sign(rand(1,1)-0.5);
    input3=unifrnd(0,bit_length)/bit_length;    sympol3=sign(rand(1,1)-0.5);
    input4=unifrnd(0,bit_length)/bit_length;    sympol4=sign(rand(1,1)-0.5);
    input5=unifrnd(0,bit_length)/bit_length;    sympol5=sign(rand(1,1)-0.5);
    
    %% FSM
    max_state=8-1;
    state=floor(max_state/2);
    rnb1=LFSR(1,bit,input1);rnb2=LFSR(2,bit,input2);rnb3=LFSR(3,bit,input3);rnb4=LFSR(4,bit,input4);rnb5=LFSR(5,bit,input5);
    for i=1:1:bit_length
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
    input(j)=(sympol1*input1+sympol2*input2+sympol3*input3+sympol4*input4+sympol5*input5);
    result(j)=output/bit_length;
    sigmoid(j)=sigmf(input(j),[12,0]);
    j
end
figure;
plot(input,result,'ro','MarkerSize',4);

hold on;
plot(input,sigmoid,'go','MarkerSize',4);
toc