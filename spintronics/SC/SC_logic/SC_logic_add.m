clear;
tic;
% input1=round(unifrnd(1,1024));
% input2=round(unifrnd(1,1024));

input1=64;
input2=256;
t=0;
for input1=0:8:256
    for input2=0:8:256
        t=t+1;
        bit_length=256;
        for j=1:1:10
            %% generate bit stream
            for i=1:1:bit_length
                rnb1(i)=s_bit(input1,bit_length);
                rnb2(i)=s_bit(input2,bit_length);
                
                if rnb1(i)==1
                    I=341;               %342
                else
                    I=85;               %86
                end
                
                if rnb2(i)==1
                    V=1;
                else
                    V=0;
                end
                out(i)=sng_function(I,V);
            end
            
            %% count
            output=sum(out);
            std_output=(input1+input2)/2;
            error(j)=output-std_output;
            error_r(j)=(output-std_output)/std_output;
        end
        sigma=mean(abs(error));
        sigma_r=mean(abs(error_r));
        in1(t)=input1;
        in2(t)=input2;
        error_mean(t)=sigma;
        error_mean_r(t)=sigma_r;
        t
    end
end
figure;
toc