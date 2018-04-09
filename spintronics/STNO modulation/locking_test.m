clear;
%% 电流关系测试 确定VCMA电压 改变电流大小
% for i=500:10:3000
%     [temp,temp2,temp3,temp4]=STNO_locking(i,0,0.5); 
%     fre(round(i/10-49))=temp(1);              %仿真频率
%     n(round(i/10-49))=i;                      %电流大小
%     mz_mean(round(i/10-49))=temp3;            %仿真mz
%     fre_predict(round(i/10-49))=temp2;        %计算的频率
%     mz_cal(round(i/10-49))=temp4;             %计算的mz
%     i
% end
%% 电压关系测试
for i=0:0.01:2.5
    [temp,temp2,temp3,temp4]=STNO_locking(2400,i,1.2 );
    fre(round(i*100+1))=temp(1);
    n(round(i*100+1))=i;
    mz_mean(round(i*100+1))=temp3;
    fre_predict(round(i*100+1))=temp2;
    mz_cal(round(i*100+1))=temp4;
    i
end

%% 画图
figure;
subplot(2,1,1);
plot(n,fre,'r');
hold on;
plot(n,fre_predict);
subplot(2,1,2);
plot(n,mz_mean,'r');
hold on;
plot(n,mz_cal);