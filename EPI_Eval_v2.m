function [] = EPI_Eval_v2(s_all,Dates,RESULTS_dir)
% get Shewart Evaluation for all time
% make text file
% make all time plots
% color based on Shewart Evaluation
% v1 showed all at once on 3by3
% v2 seperates


EPI_binary = strcmp({s_all.type}, 'EPI STABILITY');
EPIIndex = find(EPI_binary);

clear GA Dist t1 t2 HCSR SliceThick SlicePos PIU PSG LCOD SNR
clear GA Dist t1 t2 PIU PSG SNR SFNR PF Drift RDC

ind=1;
for i=1:length(EPIIndex)
    index = EPIIndex(i);
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

[t1_sorted,t1_order]=sort(t1);
Shewart_all = [];

ms=8; fs = 12; k =15;
%set(gcf, 'units','normalized','outerposition',[0.05 0.05 .9 .9]); 

%% GA
figure; 
Param = GA; rownum = 1;
Val_sorted = Param(t1_order);
Shewarts_sorted = zeros(size(Val_sorted));
for i2 = 30:length(Val_sorted)
    Val_short = Val_sorted(1:i2);
%     Shewarts_short = Shewarts_sorted(1:i2);
%     Val_short(Shewarts_short>0) = [];
    [Shewart,~] = ShewartEval(Val_short);
    Shewarts_sorted(i2) = Shewart;
end
Shewart_all(rownum,:) = Shewarts_sorted;

%subplot(3,3,1);
plot(t1_sorted,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
line([min(t1)-1 max(t1)+1], [190 190],'color',[0 0 0],'linestyle','--');
plot(t1_sorted(Shewarts_sorted==1),Val_sorted(Shewarts_sorted==1),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',ms); hold on;
plot(t1_sorted(Shewarts_sorted>1),Val_sorted(Shewarts_sorted>1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms); hold on;
ylabel('Mean Phantom Diameter'); title('Geometric Accuracy');
xlim([min(t1)-1 max(t1)+1]); 
ylim([min(0.98*min(GA),0.98*188) max(1.02*max(GA),1.02*192)])
set(gca,'FontSize',fs); grid off;

drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_All_1.png']); pause(1);
xlim([max(t1)-99 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_100_1.png']); pause(1);
xlim([max(t1)-30 max(t1)]);
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_30_1.png']); pause(1);
xlim([max(t1)-6 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_7_1.png']); pause(1);
close all

%% Dist
figure;
Param = Dist; rownum = 2;
Val_sorted = Param(t1_order);
Shewarts_sorted = zeros(size(Val_sorted));
for i2 = 30:length(Val_sorted)
    Val_short = Val_sorted(1:i2);
%     Shewarts_short = Shewarts_sorted(1:i2);
%     Val_short(Shewarts_short>0) = [];
    [Shewart,~] = ShewartEval(Val_short);
    Shewarts_sorted(i2) = Shewart;
end
Shewart_all(rownum,:) = Shewarts_sorted;

if length(Dist_sorted)<(2*k+1)
    k = floor((length(Dist_sorted)-1)/2);
end

%subplot(3,3,2);
plot(t1_sorted,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
s = movingstd(Val_sorted,k,[]); s(1:7)=s(8); s(end-7:end)=s(end-7);
M = conv(Val_sorted, ones(k,1)/k, 'same'); M(1:7)=M(8); M(end-7:end)=M(end-7);
plot(t1_sorted,M,'k-'); plot(t1_sorted,M-2*s,'k--'); plot(t1_sorted,M+2*s,'k--');
plot(t1_sorted(Shewarts_sorted==1),Val_sorted(Shewarts_sorted==1),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',ms); hold on;
plot(t1_sorted(Shewarts_sorted>1),Val_sorted(Shewarts_sorted>1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms); hold on;
ylabel('% Distortion'); title('Geometric Distortion');
xlim([min(t1)-1 max(t1)+1]); 
ylim([0 max(1.1*max(Dist),5)]);
set(gca,'FontSize',fs); grid off;

drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_All_2.png']); pause(1);
xlim([max(t1)-99 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_100_2.png']); pause(1);
xlim([max(t1)-30 max(t1)]);
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_30_2.png']); pause(1);
xlim([max(t1)-6 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_7_2.png']); pause(1);
close all
 
%% PIU
figure;
Param = PIU; 
Val_sorted = Param(t1_order);
PIU_fail = Val_sorted<.82;

%subplot(3,3,3);
plot(t1_sorted,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
line([min(t1)-1 max(t1)+1], [.82 .82],'color',[0 0 0],'linestyle','--');
plot(t1_sorted(PIU_fail==1),Val_sorted(PIU_fail==1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms); hold on;
ylabel('% Image Uniformity'); title('Image Intensity Uniformity');
xlim([min(t1)-1 max(t1)+1]); 
ylim([0 1]);
set(gca,'FontSize',fs); grid off;

drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_All_3.png']); pause(1);
xlim([max(t1)-99 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_100_3.png']); pause(1);
xlim([max(t1)-30 max(t1)]);
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_30_3.png']); pause(1);
xlim([max(t1)-6 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_7_3.png']); pause(1);
close all

%% PSG
figure;
Param = PSG; 
Val_sorted = Param(t1_order);
PSG_fail = Val_sorted>.025;

%subplot(3,3,4);
plot(t1_sorted,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
line([min(t1)-1 max(t1)+1], [.025 .025],'color',[0 0 0],'linestyle','--');
plot(t1_sorted(PSG_fail==1),Val_sorted(PSG_fail==1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms); hold on;
ylabel('% Signal Ghosting'); title('Ghosting Ratio');
xlim([min(t1)-1 max(t1)+1]); 
ylim([0 max(1.02*max(PSG),1.1*.025)]);
set(gca,'FontSize',fs); grid off;

drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_All_4.png']); pause(1);
xlim([max(t1)-99 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_100_4.png']); pause(1);
xlim([max(t1)-30 max(t1)]);
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_30_4.png']); pause(1);
xlim([max(t1)-6 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_7_4.png']); pause(1);
close all

%% SNR
figure;
Param = SNR; rownum = 3;
Val_sorted = Param(t1_order);
Shewarts_sorted = zeros(size(Val_sorted));
for i2 = 30:length(Val_sorted)
    Val_short = Val_sorted(1:i2);
    [Shewart,~] = ShewartEval(Val_short);
    Shewarts_sorted(i2) = Shewart;
end
Shewart_all(rownum,:) = Shewarts_sorted;

%subplot(3,3,5);
plot(t1_sorted,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
s = movingstd(Val_sorted,k,[]); s(1:7)=s(8); s(end-7:end)=s(end-7);
M = conv(Val_sorted, ones(k,1)/k, 'same'); M(1:7)=M(8); M(end-7:end)=M(end-7);
plot(t1_sorted,M,'k-'); plot(t1_sorted,M-2*s,'k--'); plot(t1_sorted,M+2*s,'k--');
plot(t1_sorted(Shewarts_sorted==1),Val_sorted(Shewarts_sorted==1),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',ms); hold on;
plot(t1_sorted(Shewarts_sorted>1),Val_sorted(Shewarts_sorted>1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms); hold on;
ylabel('SNR'); title('Signal To Noise Ratio');
xlim([min(t1)-1 max(t1)+1]);
set(gca,'FontSize',fs); grid off;

drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_All_5.png']); pause(1);
xlim([max(t1)-99 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_100_5.png']); pause(1);
xlim([max(t1)-30 max(t1)]);
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_30_5.png']); pause(1);
xlim([max(t1)-6 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_7_5.png']); pause(1);
close all

%% SFNR
figure;
Param = SFNR;
Val_sorted = Param(t1_order);

t1_sorted_SFNR = t1_sorted; 
t1_sorted_SFNR(isnan(Val_sorted))=[]; Val_sorted(isnan(Val_sorted))=[];
t1_sorted_SFNR(isinf(Val_sorted))=[]; Val_sorted(isinf(Val_sorted))=[];

Shewarts_sorted = zeros(size(Val_sorted));
for i2 = 30:length(Val_sorted)
    Val_short = Val_sorted(1:i2);
    [Shewart,~] = ShewartEval(Val_short);
    Shewarts_sorted(i2) = Shewart;
end
Shewart_SFNR = Shewarts_sorted;

%subplot(3,3,6);
plot(t1_sorted_SFNR,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
s = movingstd(Val_sorted,k,[]); s(1:7)=s(8); s(end-7:end)=s(end-7);
M = conv(Val_sorted, ones(k,1)/k, 'same'); M(1:7)=M(8); M(end-7:end)=M(end-7);
plot(t1_sorted_SFNR,M,'k-'); plot(t1_sorted_SFNR,M-2*s,'k--'); plot(t1_sorted_SFNR,M+2*s,'k--');
plot(t1_sorted_SFNR(Shewart_SFNR==1),Val_sorted(Shewart_SFNR==1),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',ms); hold on;
plot(t1_sorted_SFNR(Shewart_SFNR>1),Val_sorted(Shewart_SFNR>1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms); hold on;
ylabel('SFNR'); title('Signal To Flunctuation Noise Ratio');
xlim([min(t1)-1 max(t1)+1]); 
set(gca,'FontSize',fs); grid off;

drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_All_6.png']); pause(1);
xlim([max(t1)-99 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_100_6.png']); pause(1);
xlim([max(t1)-30 max(t1)]);
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_30_6.png']); pause(1);
xlim([max(t1)-6 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_7_6.png']); pause(1);
close all

%% PF
figure;
Param = PF; 
Val_sorted = Param(t1_order);

t1_sorted_PF = t1_sorted; 
t1_sorted_PF(isnan(Val_sorted))=[]; Val_sorted(isnan(Val_sorted))=[];
t1_sorted_PF(isinf(Val_sorted))=[]; Val_sorted(isinf(Val_sorted))=[];

Shewarts_sorted = zeros(size(Val_sorted));
for i2 = 30:length(Val_sorted)
    Val_short = Val_sorted(1:i2);
    [Shewart,~] = ShewartEval(Val_short);
    Shewarts_sorted(i2) = Shewart;
end
Shewart_PF = Shewarts_sorted;

%subplot(3,3,7);
plot(t1_sorted_PF,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
s = movingstd(Val_sorted,k,[]); s(1:7)=s(8); s(end-7:end)=s(end-7);
M = conv(Val_sorted, ones(k,1)/k, 'same'); M(1:7)=M(8); M(end-7:end)=M(end-7);
plot(t1_sorted_PF,M,'k-'); plot(t1_sorted_PF,M-2*s,'k--'); plot(t1_sorted_PF,M+2*s,'k--');
plot(t1_sorted_PF(Shewart_PF==1),Val_sorted(Shewart_PF==1),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',ms); hold on;
plot(t1_sorted_PF(Shewart_PF>1),Val_sorted(Shewart_PF>1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms); hold on;
ylabel('% Fluctuation'); title('Percent Fluctuation');
xlim([min(t1)-1 max(t1)+1]); 
set(gca,'FontSize',fs); grid off;
ylim([0 .4])

drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_All_7.png']); pause(1);
xlim([max(t1)-99 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_100_7.png']); pause(1);
xlim([max(t1)-30 max(t1)]);
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_30_7.png']); pause(1);
xlim([max(t1)-6 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_7_7.png']); pause(1);
close all

%% Drift
figure;
Param = Drift; 
Val_sorted = Param(t1_order);

t1_sorted_Drift = t1_sorted; 
t1_sorted_Drift(isnan(Val_sorted))=[]; Val_sorted(isnan(Val_sorted))=[];
t1_sorted_Drift(isinf(Val_sorted))=[]; Val_sorted(isinf(Val_sorted))=[];

Shewarts_sorted = zeros(size(Val_sorted));
for i2 = 30:length(Val_sorted)
    Val_short = Val_sorted(1:i2);
    [Shewart,~] = ShewartEval(Val_short);
    Shewarts_sorted(i2) = Shewart;
end
Shewart_Drift = Shewarts_sorted;

%subplot(3,3,8);
plot(t1_sorted_Drift,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
s = movingstd(Val_sorted,k,[]); s(1:7)=s(8); s(end-7:end)=s(end-7);
M = conv(Val_sorted, ones(k,1)/k, 'same'); M(1:7)=M(8); M(end-7:end)=M(end-7);
plot(t1_sorted_Drift,M,'k-'); plot(t1_sorted_Drift,M-2*s,'k--'); plot(t1_sorted_Drift,M+2*s,'k--');
plot(t1_sorted_Drift(Shewart_Drift==1),Val_sorted(Shewart_Drift==1),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',ms); hold on;
plot(t1_sorted_Drift(Shewart_Drift>1),Val_sorted(Shewart_Drift>1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms); hold on;
ylabel('Drift'); title('Drift');
xlim([min(t1)-1 max(t1)+1]); 
set(gca,'FontSize',fs); grid off;
ylim([0 2])

drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_All_8.png']); pause(1);
xlim([max(t1)-99 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_100_8.png']); pause(1);
xlim([max(t1)-30 max(t1)]);
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_30_8.png']); pause(1);
xlim([max(t1)-6 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_7_8.png']); pause(1);
close all

%% RDC
figure;
Param = RDC; 
Val_sorted = Param(t1_order);

t1_sorted_RDC = t1_sorted; 
t1_sorted_RDC(isnan(Val_sorted))=[]; Val_sorted(isnan(Val_sorted))=[];
t1_sorted_RDC(isinf(Val_sorted))=[]; Val_sorted(isinf(Val_sorted))=[];

Shewarts_sorted = zeros(size(Val_sorted));
for i2 = 30:length(Val_sorted)
    Val_short = Val_sorted(1:i2);
    [Shewart,~] = ShewartEval(Val_short);
    Shewarts_sorted(i2) = Shewart;
end
Shewart_RDC = Shewarts_sorted;

%subplot(3,3,9);
plot(t1_sorted_RDC,Val_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
s = movingstd(Val_sorted,k,[]); s(1:7)=s(8); s(end-7:end)=s(end-7);
M = conv(Val_sorted, ones(k,1)/k, 'same'); M(1:7)=M(8); M(end-7:end)=M(end-7);
plot(t1_sorted_RDC,M,'k-'); plot(t1_sorted_RDC,M-2*s,'k--'); plot(t1_sorted_RDC,M+2*s,'k--');
plot(t1_sorted_RDC(Shewart_RDC==1),Val_sorted(Shewart_RDC==1),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',ms); hold on;
plot(t1_sorted_RDC(Shewart_RDC>1),Val_sorted(Shewart_RDC>1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',ms); hold on;
ylabel('RDC'); title('Radius of Decorrelation');
xlim([min(t1)-1 max(t1)+1]); 
set(gca,'FontSize',fs); grid off;

drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_All_9.png']); pause(1);
xlim([max(t1)-99 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_100_9.png']); pause(1);
xlim([max(t1)-30 max(t1)]);
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_30_9.png']); pause(1);
xlim([max(t1)-6 max(t1)]); 
drawnow; pause(1);
print_current_figure(200,[RESULTS_dir filesep 'EPI_7_9.png']); pause(1);
close all

%% Save Figures
% All Time, 100 days, 30 days, 7 days

%% text output
% All Time, 100 days, 30 days, 7 days

fileID = fopen([RESULTS_dir filesep 'EPI_AllTime.txt'],'w');

fprintf(fileID,'Geometric Accuracy:\n');
for i = 1:length(Shewart_all(1,:))
    if Shewart_all(1,i)>0
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(1,i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Geometric Distortion:\n');
for i = 1:length(Shewart_all(2,:))
    if Shewart_all(2,i)>0
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(2,i)) '\n']);
    end
end
fprintf(fileID,'\n');
    
fprintf(fileID,'Percent Intensity Uniformity < 0.82:\n');
for i = 1:length(PIU_fail)
    if PIU_fail(i)>0
        fprintf(fileID,[datestr(t1_sorted(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Signal Ghosting > 0.025:\n');
for i = 1:length(PSG_fail)
    if PSG_fail(i)>0
        fprintf(fileID,[datestr(t1_sorted(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'SNR:\n');
for i = 1:length(Shewart_all(3,:))
    if Shewart_all(3,i)>0
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(3,i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Signal To Fluctuation Noise Ratio:\n');
for i = 1:length(Shewart_SFNR)
    if Shewart_SFNR(i)>0
        fprintf(fileID,[datestr(t1_sorted_SFNR(i)) '\tShewart Number:' num2str(Shewart_SFNR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Fluctuation:\n');
for i = 1:length(Shewart_PF)
    if Shewart_PF(i)>0
        fprintf(fileID,[datestr(t1_sorted_PF(i)) '\tShewart Number:' num2str(Shewart_PF(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Drift:\n');
for i = 1:length(Shewart_Drift)
    if Shewart_Drift(i)>0
        fprintf(fileID,[datestr(t1_sorted_Drift(i)) '\tShewart Number:' num2str(Shewart_Drift(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Radius Decorrelation:\n');
for i = 1:length(Shewart_RDC)
    if Shewart_RDC(i)>0
        fprintf(fileID,[datestr(t1_sorted_RDC(i)) '\tShewart Number:' num2str(Shewart_RDC(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fclose(fileID);

%%%%%% 100 days %%%%%%%%%

fileID = fopen([RESULTS_dir filesep 'EPI_100Days.txt'],'w');

fprintf(fileID,'Geometric Accuracy:\n');
for i = 1:length(Shewart_all(1,:))
    if Shewart_all(1,i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<100
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(1,i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Geometric Distortion:\n');
for i = 1:length(Shewart_all(2,:))
    if Shewart_all(2,i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<100
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(2,i)) '\n']);
    end
end
fprintf(fileID,'\n');
    
fprintf(fileID,'Percent Intensity Uniformity < 0.82:\n');
for i = 1:length(PIU_fail)
    if PIU_fail(i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<100
        fprintf(fileID,[datestr(t1_sorted(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Signal Ghosting > 0.025:\n');
for i = 1:length(PSG_fail)
    if PSG_fail(i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<100
        fprintf(fileID,[datestr(t1_sorted(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'SNR:\n');
for i = 1:length(Shewart_all(3,:))
    if Shewart_all(3,i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<100
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(3,i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Signal To Fluctuation Noise Ratio:\n');
for i = 1:length(Shewart_SFNR)
    if Shewart_SFNR(i)>0 && (datenum(max(t1_sorted_SFNR))-datenum(t1_sorted_SFNR(i)))<100
        fprintf(fileID,[datestr(t1_sorted_SFNR(i)) '\tShewart Number:' num2str(Shewart_SFNR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Fluctuation:\n');
for i = 1:length(Shewart_PF)
    if Shewart_PF(i)>0 && (datenum(max(t1_sorted_PF))-datenum(t1_sorted_PF(i)))<100
        fprintf(fileID,[datestr(t1_sorted_PF(i)) '\tShewart Number:' num2str(Shewart_PF(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Drift:\n');
for i = 1:length(Shewart_Drift)
    if Shewart_Drift(i)>0 && (datenum(max(t1_sorted_Drift))-datenum(t1_sorted_Drift(i)))<100
        fprintf(fileID,[datestr(t1_sorted_Drift(i)) '\tShewart Number:' num2str(Shewart_Drift(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Radius Decorrelation:\n');
for i = 1:length(Shewart_RDC)
    if Shewart_RDC(i)>0 && (datenum(max(t1_sorted_RDC))-datenum(t1_sorted_RDC(i)))<100
        fprintf(fileID,[datestr(t1_sorted_RDC(i)) '\tShewart Number:' num2str(Shewart_RDC(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fclose(fileID);

%%%%%% 30 days %%%%%%%%%

fileID = fopen([RESULTS_dir filesep 'EPI_30Days.txt'],'w');

fprintf(fileID,'Geometric Accuracy:\n');
for i = 1:length(Shewart_all(1,:))
    if Shewart_all(1,i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<30
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(1,i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Geometric Distortion:\n');
for i = 1:length(Shewart_all(2,:))
    if Shewart_all(2,i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<30
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(2,i)) '\n']);
    end
end
fprintf(fileID,'\n');
    
fprintf(fileID,'Percent Intensity Uniformity < 0.82:\n');
for i = 1:length(PIU_fail)
    if PIU_fail(i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<30
        fprintf(fileID,[datestr(t1_sorted(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Signal Ghosting > 0.025:\n');
for i = 1:length(PSG_fail)
    if PSG_fail(i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<30
        fprintf(fileID,[datestr(t1_sorted(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'SNR:\n');
for i = 1:length(Shewart_all(3,:))
    if Shewart_all(3,i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<30
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(3,i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Signal To Fluctuation Noise Ratio:\n');
for i = 1:length(Shewart_SFNR)
    if Shewart_SFNR(i)>0 && (datenum(max(t1_sorted_SFNR))-datenum(t1_sorted_SFNR(i)))<30
        fprintf(fileID,[datestr(t1_sorted_SFNR(i)) '\tShewart Number:' num2str(Shewart_SFNR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Fluctuation:\n');
for i = 1:length(Shewart_PF)
    if Shewart_PF(i)>0 && (datenum(max(t1_sorted_PF))-datenum(t1_sorted_PF(i)))<30
        fprintf(fileID,[datestr(t1_sorted_PF(i)) '\tShewart Number:' num2str(Shewart_PF(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Drift:\n');
for i = 1:length(Shewart_Drift)
    if Shewart_Drift(i)>0 && (datenum(max(t1_sorted_Drift))-datenum(t1_sorted_Drift(i)))<30
        fprintf(fileID,[datestr(t1_sorted_Drift(i)) '\tShewart Number:' num2str(Shewart_Drift(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Radius Decorrelation:\n');
for i = 1:length(Shewart_RDC)
    if Shewart_RDC(i)>0 && (datenum(max(t1_sorted_RDC))-datenum(t1_sorted_RDC(i)))<30
        fprintf(fileID,[datestr(t1_sorted_RDC(i)) '\tShewart Number:' num2str(Shewart_RDC(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fclose(fileID);

%%%%% 7 %%%%%

fileID = fopen([RESULTS_dir filesep 'EPI_7Days.txt'],'w');

fprintf(fileID,'Geometric Accuracy:\n');
for i = 1:length(Shewart_all(1,:))
    if Shewart_all(1,i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<7
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(1,i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Geometric Distortion:\n');
for i = 1:length(Shewart_all(2,:))
    if Shewart_all(2,i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<7
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(2,i)) '\n']);
    end
end
fprintf(fileID,'\n');
    
fprintf(fileID,'Percent Intensity Uniformity < 0.82:\n');
for i = 1:length(PIU_fail)
    if PIU_fail(i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<7
        fprintf(fileID,[datestr(t1_sorted(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Signal Ghosting > 0.025:\n');
for i = 1:length(PSG_fail)
    if PSG_fail(i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<7
        fprintf(fileID,[datestr(t1_sorted(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'SNR:\n');
for i = 1:length(Shewart_all(3,:))
    if Shewart_all(3,i)>0 && (datenum(max(t1_sorted))-datenum(t1_sorted(i)))<7
        fprintf(fileID,[datestr(t1_sorted(i)) '\tShewart Number:' num2str(Shewart_all(3,i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Signal To Fluctuation Noise Ratio:\n');
for i = 1:length(Shewart_SFNR)
    if Shewart_SFNR(i)>0 && (datenum(max(t1_sorted_SFNR))-datenum(t1_sorted_SFNR(i)))<7
        fprintf(fileID,[datestr(t1_sorted_SFNR(i)) '\tShewart Number:' num2str(Shewart_SFNR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Fluctuation:\n');
for i = 1:length(Shewart_PF)
    if Shewart_PF(i)>0 && (datenum(max(t1_sorted_PF))-datenum(t1_sorted_PF(i)))<7
        fprintf(fileID,[datestr(t1_sorted_PF(i)) '\tShewart Number:' num2str(Shewart_PF(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Percent Drift:\n');
for i = 1:length(Shewart_Drift)
    if Shewart_Drift(i)>0 && (datenum(max(t1_sorted_Drift))-datenum(t1_sorted_Drift(i)))<7
        fprintf(fileID,[datestr(t1_sorted_Drift(i)) '\tShewart Number:' num2str(Shewart_Drift(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Radius Decorrelation:\n');
for i = 1:length(Shewart_RDC)
    if Shewart_RDC(i)>0 && (datenum(max(t1_sorted_RDC))-datenum(t1_sorted_RDC(i)))<7
        fprintf(fileID,[datestr(t1_sorted_RDC(i)) '\tShewart Number:' num2str(Shewart_RDC(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fclose(fileID);


