clear;
iris_data = csvread('iris_test.csv',1,0);
layer1_weights=[6.9032,	5.04648,	-4.34227,	-5.02293;
            -6.24207,	-5.23362,	5.34296,	3.99471;
            14.8788,	13.9252,	-13.9821,	-14.523;
            11.3529,	12.4055,	-11.7516,	-12.2361];

layer1_biases=[-0.889334,   -0.844831,  0.799981,   0.926365];

logits_weights=[-13.7329,	1.27417,	14.1756;
            -13.4678,	-1.06278,	13.8347;
            12.981,	1.1352,	-13.7875;
            12.9109,	1.89396,	-13.8149];

logits_biases=[0.0456089,   0.852301,   -0.779679];

layer1_weights=layer1_weights./10;
logits_weights=logits_weights./15;
layer1_biases=layer1_biases./1;
logits_biases=logits_biases./15;

accuracy=0;
for i=1:1:30
%     i=24;
    test=0.1*iris_data(i,1:4);
    layer1_pre=test*layer1_weights+layer1_biases;
    layer1=sigmf(layer1_pre,[1,0]);
    output=layer1*logits_weights+logits_biases;
    [temp,index]=max(output);
    if (index==iris_data(i,5)+1)
        accuracy=accuracy+1;
    end
end