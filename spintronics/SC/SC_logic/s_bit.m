function [ bit ] = s_bit( input  , bits_length )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    random_number=round(unifrnd(0,bits_length));
    input=round(input);
    if input>random_number
        bit=1;
    else
        bit=0;
    end

end

