clear;
iris_data = csvread('iris_test.csv',1,0);
% layer1_weights=[6.9032,	5.04648,	-4.34227,	-5.02293;
%             -6.24207,	-5.23362,	5.34296,	3.99471;
%             14.8788,	13.9252,	-13.9821,	-14.523;
%             11.3529,	12.4055,	-11.7516,	-12.2361];
%
% layer1_biases=[-0.889334,   -0.844831,  0.799981,   0.926365];
%
% logits_weights=[-13.7329,	1.27417,	14.1756;
%             -13.4678,	-1.06278,	13.8347;
%             12.981,	1.1352,	-13.7875;
%             12.9109,	1.89396,	-13.8149];
%
% logits_biases=[0.0456089,   0.852301,   -0.779679];


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

layer1_weights=layer1_weights./12;
logits_weights=logits_weights./12;
layer1_biases=layer1_biases./12;
logits_biases=logits_biases./12;

accuracy=0;
for i=1:1:30
    %     i=24;
    test=0.1*iris_data(i,1:4);
    layer1_pre=test*layer1_weights+layer1_biases;
    layer1=sigmf(layer1_pre,[12,0]);
    output=layer1*logits_weights+logits_biases;
    [temp,index]=max(output);
    if (index==iris_data(i,5)+1)
        accuracy=accuracy+1;
    end
end