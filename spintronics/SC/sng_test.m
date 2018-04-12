% clear;
tic;
% for I=560:1:560
%     positive=0;
%     negtive=0;
%     for i=1:1:8192*2
%         mz=sng_function(I,0);
%         if(mz==1)
%             positive=positive+1;
%         else
%             negtive=negtive+1;
%         end
%     end
%     I
%     p(I)=positive/(positive+negtive);
%     p(I)
% end
% figure;
% plot(p);

% k=1;
% for V=-1:0.02:1
%     positive=0;
%     negtive=0;
%     for i=1:1:4096
%         mz=sng_function(340,V);
%         if(mz>0)
%             positive=positive+1;
%         else
%             negtive=negtive+1;
%         end
%     end
%     V
%     Voltage(k)=V;
%     p(k)=positive/(positive+negtive);
%     p(k)
%     k=k+1;
% end
% figure;
% plot(Voltage,p);
% toc

 for i=-1:0.02:1
    y(round(50*i+51))=sigmf(i,[6,0.03]);
 end
figure;
plot(Voltage,p);
hold on
plot(Voltage,y)