clear;
tic;
input1=round(unifrnd(1,256));
input2=round(unifrnd(1,256));

t=0;
for input1=0:16:256
    for input2=0:16:256
        t=t+1;
        bit_length=256;
        for j=1:1:10
            %% generate bit stream
            for i=1:1:bit_length
                rnb1(i)=s_bit(input1,bit_length);
                rnb2(i)=s_bit(input2,bit_length);
                
                scl_rn(i)=round(unifrnd(0,1));
                if scl_rn(i)==1
                    out(i)=rnb1(i);
                else
                    out(i)=rnb2(i);
                end
            end
            
            %% count
            output=sum(out);
            std_output=(input1+input2)/2;
            error(j)=output-std_output;
            
        end
        sigma=mean(abs(error));
        in1(t)=input1;
        in2(t)=input2;
        error_mean(t)=sigma;
        t
    end
end
figure;
toc