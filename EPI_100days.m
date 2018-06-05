function []=EPI_100days(s_all,Dates,RESULTS_dir)

% grab all with s_all.type = 'EPI'
EPI_binary = strcmp({s_all.type}, 'EPI STABILITY');
EPIIndex = find(EPI_binary);
% get most recent date, 100 days, all time dates
EPIMostRecent = max(Dates(EPI_binary==1));
EPIMostRecentIndex = find(Dates==EPIMostRecent & EPI_binary');

EPI100Index = find(EPI_binary' .* (EPIMostRecent-Dates<100));
%T1100Index = find(Dates(T1_binary' & T1MostRecent-Dates<100));
%T1_100 = Dates(T1_index' & T1MostRecent-Dates<100);
clear GA Dist t1 t2 HCSR SliceThick SlicePos PIU PSG LCOD SNR

% make 3 by 3 subplot of values
ind=1;
for i=1:length(EPI100Index)
    index = EPI100Index(i);
    %if sum(strfind(s_all(index).name, 'ACR'))==0;
    GA(ind) = mean(s_all(index).GeometricDistortion(1:4));
    Dist(ind) = s_all(index).GeometricDistortion(5);
    t1(ind) = datetime(s_all(index).scan_date,'InputFormat','yyyyMMdd');
    t2(ind) = datenum(datetime(s_all(index).scan_date,'InputFormat','yyyyMMdd'));
    PIU(ind) = s_all(index).PIU(1);
    PSG(ind) = s_all(index).PSG(1);
    SNR(ind) = s_all(index).SNR;
    SFNR(ind) = s_all(index).SFNR;
    PF(ind) = s_all(index).PF;
    Drift(ind) = s_all(index).Drift;
    RDC(ind) = s_all(index).RDC;
    ind=ind+1;
    %end
end

figure; set(gcf, 'units','normalized','outerposition',[0.05 0.05 .9 .9]); 
ms=8; fs = 12; k =15;

subplot(3,3,1);
plot(t1,GA,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(GA)
%     if abs(190-GA(i))>2; col = [1 0 0]; else col = [0 1 0]; end
%     plot(t1(i),GA(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
% end
line([min(t2)-1 max(t2)+1], [190 190],'color',[0 0 0],'linestyle','--');
ylabel('Mean Phantom Diameter'); title('Geometric Accuracy');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(0.98*min(GA),0.98*188) max(1.02*max(GA),1.02*192)])
set(gca,'FontSize',fs); grid off;
textGA = t1(GA<188 | GA>192);

subplot(3,3,2);
[t1_sorted,t1_order]=sort(t1);
Dist_sorted = Dist(t1_order);
plot(t1_sorted,Dist_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),Dist_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
s = movingstd(Dist_sorted,k,[]);
M = conv(Dist_sorted, ones(k,1)/k, 'same');
s(1:7)=s(8);s(end-7:end)=s(end-7);
M(1:7)=M(8);M(end-7:end)=M(end-7);
plot(t1_sorted,M,'k-');
plot(t1_sorted,M-2*s,'k--');
plot(t1_sorted,M+2*s,'k--');
ylabel('% Distortion'); title('Geometric Distortion');
xlim([min(t2)-1 max(t2)+1]); 
ylim([0 max(1.1*max(Dist),5)]);
set(gca,'FontSize',fs); grid off;
textDist = t1(Dist_sorted>M+2*s);

subplot(3,3,3);
[t1_sorted,t1_order]=sort(t1);
Val_sorted = PIU(t1_order);
plot(t1_sorted,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),Dist_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
line([min(t2)-1 max(t2)+1], [.82 .82],'color',[0 0 0],'linestyle','--');
ylabel('% Image Uniformity'); title('Image Intensity Uniformity');
xlim([min(t2)-1 max(t2)+1]); 
ylim([0 1]);
set(gca,'FontSize',fs); grid off;
textPIU = t1(PIU<.82);

subplot(3,3,4);
[t1_sorted,t1_order]=sort(t1);
Val_sorted = PSG(t1_order);
plot(t1_sorted,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),Dist_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
line([min(t2)-1 max(t2)+1], [.025 .025],'color',[0 0 0],'linestyle','--');
ylabel('% Signal Ghosting'); title('Ghosting Ratio');
xlim([min(t2)-1 max(t2)+1]); 
ylim([0 max(1.02*max(PSG),1.1*.025)]);
set(gca,'FontSize',fs); grid off;
textG = t1(PSG>.025);

subplot(3,3,5);
[t1_sorted,t1_order]=sort(t1);
Val_sorted = SNR(t1_order);
plot(t1_sorted,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),Dist_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
s = movingstd(Val_sorted,k,[]);
M = conv(Val_sorted, ones(k,1)/k, 'same');
s(1:7)=s(8);s(end-7:end)=s(end-7);
M(1:7)=M(8);M(end-7:end)=M(end-7);
plot(t1_sorted,M,'k-');
plot(t1_sorted,M-2*s,'k--');
plot(t1_sorted,M+2*s,'k--');
ylabel('SNR'); title('Signal To Noise Ratio');
xlim([min(t2)-1 max(t2)+1]);
set(gca,'FontSize',fs); grid off;

subplot(3,3,6);
[t1_sorted,t1_order]=sort(t1);
Val_sorted = SFNR(t1_order);
t1_sorted_curr = t1_sorted; 
t1_sorted_curr(isnan(Val_sorted))=[]; Val_sorted(isnan(Val_sorted))=[];
t1_sorted_curr(isinf(Val_sorted))=[]; Val_sorted(isinf(Val_sorted))=[];
plot(t1_sorted_curr,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),Dist_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
s = movingstd(Val_sorted,k,[]);
M = conv(Val_sorted, ones(k,1)/k, 'same');
s(1:7)=s(8);s(end-7:end)=s(end-7);
M(1:7)=M(8);M(end-7:end)=M(end-7);
plot(t1_sorted_curr,M,'k-');
plot(t1_sorted_curr,M-2*s,'k--');
plot(t1_sorted_curr,M+2*s,'k--');
ylabel('SFNR'); title('Signal To Flunctuation Noise Ratio');
xlim([min(t2)-1 max(t2)+1]); 
set(gca,'FontSize',fs); grid off;

subplot(3,3,7);
[t1_sorted,t1_order]=sort(t1);
Val_sorted = PF(t1_order);
t1_sorted_curr = t1_sorted; 
t1_sorted_curr(isnan(Val_sorted))=[]; Val_sorted(isnan(Val_sorted))=[];
t1_sorted_curr(isinf(Val_sorted))=[]; Val_sorted(isinf(Val_sorted))=[];
plot(t1_sorted_curr,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),Dist_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
s = movingstd(Val_sorted,k,[]);
M = conv(Val_sorted, ones(k,1)/k, 'same');
s(1:7)=s(8);s(end-7:end)=s(end-7);
M(1:7)=M(8);M(end-7:end)=M(end-7);
plot(t1_sorted_curr,M,'k-');
plot(t1_sorted_curr,M-2*s,'k--');
plot(t1_sorted_curr,M+2*s,'k--');
ylabel('% Fluctuation'); title('Percent Fluctuation');
xlim([min(t2)-1 max(t2)+1]); 
set(gca,'FontSize',fs); grid off;

subplot(3,3,8);
[t1_sorted,t1_order]=sort(t1);
Val_sorted = Drift(t1_order);
t1_sorted_curr = t1_sorted; 
t1_sorted_curr(isnan(Val_sorted))=[]; Val_sorted(isnan(Val_sorted))=[];
t1_sorted_curr(isinf(Val_sorted))=[]; Val_sorted(isinf(Val_sorted))=[];
plot(t1_sorted_curr,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),Dist_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
s = movingstd(Val_sorted,k,[]);
M = conv(Val_sorted, ones(k,1)/k, 'same');
s(1:7)=s(8);s(end-7:end)=s(end-7);
M(1:7)=M(8);M(end-7:end)=M(end-7);
plot(t1_sorted_curr,M,'k-');
plot(t1_sorted_curr,M-2*s,'k--');
plot(t1_sorted_curr,M+2*s,'k--');
ylabel('Drift'); title('Drift');
xlim([min(t2)-1 max(t2)+1]); 
set(gca,'FontSize',fs); grid off;

subplot(3,3,9);
[t1_sorted,t1_order]=sort(t1);
Val_sorted = RDC(t1_order);
t1_sorted_curr = t1_sorted; 
t1_sorted_curr(isnan(Val_sorted))=[]; Val_sorted(isnan(Val_sorted))=[];
t1_sorted_curr(isinf(Val_sorted))=[]; Val_sorted(isinf(Val_sorted))=[];
plot(t1_sorted_curr,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),Dist_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
s = movingstd(Val_sorted,k,[]);
M = conv(Val_sorted, ones(k,1)/k, 'same');
s(1:7)=s(8);s(end-7:end)=s(end-7);
M(1:7)=M(8);M(end-7:end)=M(end-7);
plot(t1_sorted_curr,M,'k-');
plot(t1_sorted_curr,M-2*s,'k--');
plot(t1_sorted_curr,M+2*s,'k--');
ylabel('RDC'); title('Radius of Decorrelation');
xlim([min(t2)-1 max(t2)+1]); 
set(gca,'FontSize',fs); grid off;

pause(2);
saveas(gcf,[RESULTS_dir filesep 'LongTermEPI.fig'])
set(gcf, 'PaperPositionMode', 'auto');
print_current_figure(400,[RESULTS_dir filesep 'LongTermEPI.png'])
pause(2); close all;



