clear;
% experiment=[998,1080,1215,1280,1365,1440,1492,1522,1604,1634];
% experiment_fre=[6.18764,6.0595,5.76659,5.675,5.474,5.29,5.144,5.016,4.8513,4.6316];
% experiment2=[1061,1121,1188,1255,1322,1389,1456,1516,1598,1687];
% experiment_fre2=[5.41876,5.29,5.144,4.9977,4.88787,4.79634,4.6865,4.595,4.5,4.37529];
% for i=1200:50:2300
%     temp=oscillator_perpendicular(i,0);
%     fre(round(i/50-23))=temp(1);
%     n(round(i/50-23))=i;
%     i
% end
% figure;
% plot(n,fre,'-s');
% hold on;
% plot(experiment,experiment_fre,'-gs');
% hold on;
% plot(experiment2,experiment_fre2,'-rs');
for i=0.1:0.1:4
    temp=oscillator_perpendicular(1000,i);
    fre(round(i/0.1))=temp(1);
    n(round(i/0.1))=i;
    i
end
% for i=1:1:99
%     fre_dre(i)=fre(i+1)-fre(i);
% end
% for i=1:1:98
%     fre_dre2(i)=fre_dre(i+1)-fre_dre(i);
% end
figure;
plot(n,fre,'-s');
% figure;
% plot(n(1:80),fre_dre(1:80));
% figure;
% plot(n(1:80),fre_dre2(1:80));