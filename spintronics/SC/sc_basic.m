clear;
tic;
input1=500;
input2=500;

bit_length=1024;

%% generate bit stream
for i=1:1:bit_length
    random_number1=unifrnd(0,bit_length);
    random_number2=unifrnd(0,bit_length);
    if input1>random_number1
        rnb1(i)=1;
    else
        rnb1(i)=0;
    end
    
    if input2>random_number2
        rnb2(i)=1;
    else
        rnb2(i)=0;
    end
    
    if rnb1(i)==0
        I=342;               %342
    else
        I=86;               %86
    end
    
    if rnb2(i)==0
        V=1;
    else
        V=0;
    end
    out(i)=sng_function(I,V);
end

%% count
output=0;
for i=1:1:bit_length
    if out(i)==1
        output=output+1;
    end
    
end
toc