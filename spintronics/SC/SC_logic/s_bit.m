function [ bit ] = s_bit( input  , bits_length )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    random_number=round(unifrnd(0,bits_length));
    input=round(input);
    if input>random_number
        bit=1;
    else
        bit=0;
    end

end

