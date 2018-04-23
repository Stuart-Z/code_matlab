clear;
tic;
bit=10;
bit_length=64;
sequence_length=256*8;
for j=1:1:400
    %%
    input1=unifrnd(0,bit_length);    sympol1=sign(rand(1,1)-0.5);
    input2=unifrnd(0,bit_length);    sympol2=sign(rand(1,1)-0.5);
    input3=unifrnd(0,bit_length);    sympol3=sign(rand(1,1)-0.5);
    input4=unifrnd(0,bit_length);    sympol4=sign(rand(1,1)-0.5);
    input5=unifrnd(0,bit_length);    sympol5=sign(rand(1,1)-0.5);
    
    
    %% generate bit stream
    voltage_c=0;
    delta_v=0.2;
    for i=1:1:sequence_length
        rnb1(i)=s_bit(input1,bit_length);
        rnb2(i)=s_bit(input2,bit_length);
        rnb3(i)=s_bit(input3,bit_length);
        rnb4(i)=s_bit(input4,bit_length);
        rnb5(i)=s_bit(input5,bit_length);
        
        voltage_c=(sympol1*rnb1(i)+sympol2*rnb2(i)+sympol3*rnb3(i)+sympol4*rnb4(i)+sympol5*rnb5(i))*delta_v+voltage_c;
        if(voltage_c<-0.9) voltage_c=-0.9;  end
        if(voltage_c>0.9) voltage_c=0.9;    end
        
%         out(i)=sng_function_nn(340,voltage_c);
        out(i)=s_bit(bit_length*sigmf(voltage_c,[6,0]),bit_length);
    end

    %% count
    output=sum(out);
    input(j)=(sympol1*input1+sympol2*input2+sympol3*input3+sympol4*input4+sympol5*input5)/bit_length;
    result(j)=output/sequence_length;
    sigmoid(j)=sigmf(input(j),[10,0]);
    j
end
figure;
plot(input,result,'ro','MarkerSize',4);

hold on;
plot(input,sigmoid,'go','MarkerSize',4);
toc