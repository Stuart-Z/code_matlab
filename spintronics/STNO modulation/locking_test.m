clear;
%% ������ϵ���� ȷ��VCMA��ѹ �ı������С
% for i=500:10:3000
%     [temp,temp2,temp3,temp4]=STNO_locking(i,0,0.5); 
%     fre(round(i/10-49))=temp(1);              %����Ƶ��
%     n(round(i/10-49))=i;                      %������С
%     mz_mean(round(i/10-49))=temp3;            %����mz
%     fre_predict(round(i/10-49))=temp2;        %�����Ƶ��
%     mz_cal(round(i/10-49))=temp4;             %�����mz
%     i
% end
%% ��ѹ��ϵ����
for i=0:0.01:2.5
    [temp,temp2,temp3,temp4]=STNO_locking(2400,i,1.2 );
    fre(round(i*100+1))=temp(1);
    n(round(i*100+1))=i;
    mz_mean(round(i*100+1))=temp3;
    fre_predict(round(i*100+1))=temp2;
    mz_cal(round(i*100+1))=temp4;
    i
end

%% ��ͼ
figure;
subplot(2,1,1);
plot(n,fre,'r');
hold on;
plot(n,fre_predict);
subplot(2,1,2);
plot(n,mz_mean,'r');
hold on;
plot(n,mz_cal);