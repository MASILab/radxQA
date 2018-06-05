% open .mat of results
% make plots, place in RESULTS
% (1) single time point measures
% (2) short term stability results
% (3) long term stability results
% for T1 axial, for EPI STABILITY

%% Preliminaries

QA_dir = ['/Volumes' filesep 'GRAID' filesep 'BlazeFake' filesep 'QA'];
NIFTI_dir = [QA_dir filesep 'NIFTIS'];
OUT_dir = [QA_dir filesep 'PROCESSED'];
RESULTS_dir = [QA_dir filesep 'RESULTS'];
addpath(genpath([QA_dir filesep 'EXTRA']));

% load s_all
load([OUT_dir filesep 'Results_ALL.mat']);

% get all dates
Dates = zeros(size(s_all,2),1);
for i=1:size(s_all,2)
    Dates(i) = datenum(datetime(s_all(i).scan_date,'InputFormat','yyyyMMdd'))
    %Dates(i) = datenum(s_all(i).scan_date,'yyyyMMdd'); % this doesn't handle year change well?
end

%% T1 axial - long term

% get dates

% grab all with s_all.type = 'T1'
T1_binary = strcmp({s_all.type}, 'T1')
T1Index = find(T1_binary);

% T1: get most recent date, 100 days, all time dates
T1MostRecent = max(Dates(T1_binary==1));
T1MostRecentIndex = find(Dates==T1MostRecent & T1_binary');

T1100Index = find(T1_binary' .* (T1MostRecent-Dates<100))
%T1100Index = find(Dates(T1_binary' & T1MostRecent-Dates<100));
%T1_100 = Dates(T1_index' & T1MostRecent-Dates<100);
clear GA Dist t1 t2 HCSR SliceThick SlicePos PIU PSG LCOD SNR

% make 3 by 3 subplot of values
ind=1;
for i=1:length(T1100Index)
    index = T1100Index(i);
    if sum(strfind(s_all(index).name, 'ACR'))==0;
        GA(ind) = mean(s_all(index).GeometricDistortion(1:4));
        Dist(ind) = s_all(index).GeometricDistortion(5);
        t1(ind) = datetime(s_all(index).scan_date,'InputFormat','yyyyMMdd');
        t2(ind) = datenum(datetime(s_all(index).scan_date,'InputFormat','yyyyMMdd'));
        HCSR(ind) = max(s_all(index).HCSpatialRes);
        SliceThick(ind) = s_all(index).SliceThick(1);
        SlicePos(ind) = s_all(index).SlicePositionS1(2);
        PIU(ind) = s_all(index).PIU(1);
        PSG(ind) = s_all(index).PSG(1);
        LCOD(ind) = s_all(index).LowContrastDet;
        SNR(ind) = s_all(index).SNR;
        ind=ind+1;
    end
end

figure; set(gcf, 'units','normalized','outerposition',[0 0 1 1]); 
ms=8; fs = 12; k =15;

subplot(3,3,1);
for i=1:length(GA)
    if abs(190-GA(i))>2; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),GA(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [188 188],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [192 192],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [190 190],'color',[0 0 0],'linestyle','-');
ylabel('Mean Phantom Diameter'); title('Geometric Accuracy');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(0.98*min(GA),0.98*188) max(1.02*max(GA),1.02*192)])
set(gca,'FontSize',fs); grid off;

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

subplot(3,3,3); % worse case HSR
for i=1:length(HCSR)
    if HCSR(i)>1; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),HCSR(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [1.05 1.05],'color',[1 0 0],'linestyle','--');
ylabel('Worst-Case Resolution'); title('High Contrast Spatial Resolution');
xlim([min(t2)-1 max(t2)+1]); 
ylim([.89 1.12])
set(gca,'FontSize',fs); grid off;

subplot(3,3,4);
for i=1:length(SliceThick)
    if abs(5-SliceThick(i))>.7; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),SliceThick(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [4.3 4.3],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [5.7 5.7],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [5 5],'color',[0 0 0],'linestyle','-');
ylabel('Slice Thickness'); title('Slice Thickness Accuracy');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(0.98*min(SliceThick),0.98*4.2) max(1.02*max(SliceThick),1.02*5.8)])
set(gca,'FontSize',fs); grid off;

subplot(3,3,5);
for i=1:length(SlicePos)
    if abs(0-SlicePos(i))>2.5; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),SlicePos(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [-2.5 -2.5],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [2.5 2.5],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [0 0],'color',[0 0 0],'linestyle','-');
ylabel('Slice Position'); title('Slice Position Accuracy');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(1.1*min(SlicePos),1.1*-2.5) max(1.1*max(SlicePos),1.1*2.5)])
set(gca,'FontSize',fs); grid off;

subplot(3,3,6);
for i=1:length(PIU)
    if PIU(i)<.82; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),PIU(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [.82 .82],'color',[1 0 0],'linestyle','--');
ylabel('% Integral Uniformity'); title('Image Intensity Uniformity');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(0.98*min(PIU),0.90*.82) max(1.02*max(PIU),1.1*.82)])
set(gca,'FontSize',fs); grid off;

subplot(3,3,7);
for i=1:length(PSG)
    if PSG(i)>.025; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),PSG(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [.025 .025],'color',[1 0 0],'linestyle','--');
ylabel('% Signal Ghosting'); title('Ghosting Ratio');
xlim([min(t2)-1 max(t2)+1]); 
ylim([0 max(1.02*max(PSG),1.1*.025)])
set(gca,'FontSize',fs); grid off;

subplot(3,3,8);
for i=1:length(LCOD)
    if LCOD(i)<37; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),LCOD(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [36.5 36.5],'color',[0 0 0],'linestyle','--');
ylabel('# Spokes Detected'); title('Low Contrast Object Detectability');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(0.98*min(LCOD),0.90*37) 40])
set(gca,'FontSize',fs); grid off;

subplot(3,3,9);
[t1_sorted,t1_order]=sort(t1);
SNR_sorted = SNR(t1_order);
plot(t1_sorted,SNR_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),SNR_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
s = movingstd(SNR_sorted,k,[]);
M = conv(SNR_sorted, ones(k,1)/k, 'same');
s(1:7)=s(8);s(end-7:end)=s(end-7);
M(1:7)=M(8);M(end-7:end)=M(end-7);
plot(t1_sorted,M,'k-');
plot(t1_sorted,M-2*s,'k--');
plot(t1_sorted,M+2*s,'k--');
ylabel('SNR'); title('Signal To Noise Ratio');
xlim([min(t2)-1 max(t2)+1]); 
%ylim([0 max(1.1*max(Dist),5)]);
set(gca,'FontSize',fs); grid off;

saveas(gcf,[RESULTS_dir filesep 'LongTermT1.fig'])
print_current_figure(400,[RESULTS_dir filesep 'LongTermT1.png'])
pause(2); close all;

%% T1-axial - single time point

name = s_all(T1MostRecentIndex).name;

load([OUT_dir filesep name filesep 'Results.mat']) % load S

A = imread([OUT_dir filesep name filesep 'SNR.png']);
B = imread([OUT_dir filesep name filesep 'ContrastWorse.png']);
C = imread([OUT_dir filesep name filesep 'ContrastBest.png']);
D = imread([OUT_dir filesep name filesep 'Ghosting.png']);
E = imread([OUT_dir filesep name filesep 'GeometricDistortion.png']);
F = imread([OUT_dir filesep name filesep 'GeometricAccuracy.png']);

if (sum(s.GeometricAccuracy<188) + sum(s.GeometricAccuracy>192)) > 0; pfg = 'FAIL'; else pfg = 'PASS'; end   
if (max(s.HCSpatialRes)) > 1; pfhc = 'FAIL'; else pfhc = 'PASS'; end   
if s.SliceThick(2) ==0; pfst = 'FAIL'; else pfst = 'PASS'; end      
if s.SlicePositionS1(3) ==0; pfsp = 'FAIL'; else pfsp = 'PASS'; end      
if s.PIU(1) < .82; pfpiu = 'FAIL'; else pfpiu = 'PASS'; end      
if s.PSG(2) ==0; pfpg = 'FAIL'; else pfpg = 'PASS'; end      
if s.LowContrastDet < 37; pflc = 'FAIL'; else pflc = 'PASS'; end      
[aa bb] =max(abs(s.GeometricAccuracy-190)); 

data = {'Scan Name',s.name,'N/A'; ...
    'Scan Type',s.type ,'N/A'; ...
    'Scanner',s.scanner, 'N/A'; ...
    'Date',s.scan_date, 'N/A'; ...
    'Description',s.series_description, 'N/A'; ...
    'Resolution (xy)',num2str(s.reso(1)), 'N/A'; ...
    'Geometric Accuracy',num2str(s.GeometricAccuracy(bb)), pfg; ...
    'Geometric Distrotion',num2str(s.GeometricDistortion(5)), 'N/A'; ...
    'Spatial Res',num2str(max(s.HCSpatialRes)), pfhc; ...
    'Slice Thickness',num2str(s.SliceThick(1)), pfst; ...
    'Slice Position',num2str(s.SlicePositionS1(2)), pfsp; ...
    '% Uniformity',num2str(s.PIU(1)), pfpiu; ...
    '% Signal Ghosting',num2str(s.PSG(1)), pfpg; ...
    '# spokes detected',num2str(s.LowContrastDet), pflc; ...
    'SNR',num2str(s.SNR), 'N/A'}   ;   

close all;
f = figure(5); set(gcf, 'units','normalized','outerposition',[0 0 1 1]); 
t = uitable(f,'Data',data,'ColumnWidth',{100,300,70})
subplot(1,2,2),plot(3)
pos = get(subplot(1,2,2),'position');
delete(subplot(2,1,2))
set(t,'units','normalized')
set(t,'position',pos)
axis off; 
subplot(3,4,1); imagesc(A); axis equal; axis off; title('SNR ROIs');
subplot(3,4,2); imagesc(D); axis equal; axis off; title('Ghosting ROIs');
subplot(3,4,5); imagesc(B); axis equal; axis off; colormap gray; title('Worst Case Contrast');
subplot(3,4,6); imagesc(C); axis equal; axis off; colormap gray; title('Best Case Contrast');
subplot(3,4,9); imagesc(E); axis equal; axis off; title('Distortion Measures');
subplot(3,4,10); imagesc(F); axis equal; axis off; title('Geometric Accuracy');

saveas(gcf,[RESULTS_dir filesep 'SingleTimeT1.fig'])
print_current_figure(400,[RESULTS_dir filesep 'SingleTimeT1.png'])
pause(2); close all;

%% T1 axial - all time

clear GA Dist t1 t2 HCSR SliceThick SlicePos PIU PSG LCOD SNR

ind=1;
for i=1:length(T1Index)
    index = T1Index(i);
    if sum(strfind(s_all(index).name, 'ACR'))==0;
        GA(ind) = mean(s_all(index).GeometricDistortion(1:4));
        Dist(ind) = s_all(index).GeometricDistortion(5);
        t1(ind) = datetime(s_all(index).scan_date,'InputFormat','yyyyMMdd');
        t2(ind) = datenum(datetime(s_all(index).scan_date,'InputFormat','yyyyMMdd'));
        HCSR(ind) = max(s_all(index).HCSpatialRes);
        SliceThick(ind) = s_all(index).SliceThick(1);
        SlicePos(ind) = s_all(index).SlicePositionS1(2);
        PIU(ind) = s_all(index).PIU(1);
        PSG(ind) = s_all(index).PSG(1);
        LCOD(ind) = s_all(index).LowContrastDet;
        SNR(ind) = s_all(index).SNR;
        ind=ind+1;
    end
end

figure; set(gcf, 'units','normalized','outerposition',[0 0 1 1]); 
ms=8; fs = 12; k =15;

subplot(3,3,1);
for i=1:length(GA)
    if abs(190-GA(i))>2; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),GA(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [188 188],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [192 192],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [190 190],'color',[0 0 0],'linestyle','-');
ylabel('Mean Phantom Diameter'); title('Geometric Accuracy');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(0.98*min(GA),0.98*188) max(1.02*max(GA),1.02*192)])
set(gca,'FontSize',fs); grid off;

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

subplot(3,3,3); % worse case HSR
for i=1:length(HCSR)
    if HCSR(i)>1; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),HCSR(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [1.05 1.05],'color',[1 0 0],'linestyle','--');
ylabel('Worst-Case Resolution'); title('High Contrast Spatial Resolution');
xlim([min(t2)-1 max(t2)+1]); 
ylim([.89 1.12])
set(gca,'FontSize',fs); grid off;

subplot(3,3,4);
for i=1:length(SliceThick)
    if abs(5-SliceThick(i))>.7; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),SliceThick(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [4.3 4.3],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [5.7 5.7],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [5 5],'color',[0 0 0],'linestyle','-');
ylabel('Slice Thickness'); title('Slice Thickness Accuracy');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(0.98*min(SliceThick),0.98*4.2) max(1.02*max(SliceThick),1.02*5.8)])
set(gca,'FontSize',fs); grid off;

subplot(3,3,5);
for i=1:length(SlicePos)
    if abs(0-SlicePos(i))>2.5; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),SlicePos(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [-2.5 -2.5],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [2.5 2.5],'color',[0 0 0],'linestyle','--');
line([min(t2)-1 max(t2)+1], [0 0],'color',[0 0 0],'linestyle','-');
ylabel('Slice Position'); title('Slice Position Accuracy');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(1.1*min(SlicePos),1.1*-2.5) max(1.1*max(SlicePos),1.1*2.5)])
set(gca,'FontSize',fs); grid off;

subplot(3,3,6);
for i=1:length(PIU)
    if PIU(i)<.82; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),PIU(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [.82 .82],'color',[1 0 0],'linestyle','--');
ylabel('% Integral Uniformity'); title('Image Intensity Uniformity');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(0.98*min(PIU),0.90*.82) max(1.02*max(PIU),1.1*.82)])
set(gca,'FontSize',fs); grid off;

subplot(3,3,7);
for i=1:length(PSG)
    if PSG(i)>.025; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),PSG(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [.025 .025],'color',[1 0 0],'linestyle','--');
ylabel('% Signal Ghosting'); title('Ghosting Ratio');
xlim([min(t2)-1 max(t2)+1]); 
ylim([0 max(1.02*max(PSG),1.1*.025)])
set(gca,'FontSize',fs); grid off;

subplot(3,3,8);
for i=1:length(LCOD)
    if LCOD(i)<37; col = [1 0 0]; else col = [0 1 0]; end
    plot(t1(i),LCOD(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',ms); hold on;
end
line([min(t2)-1 max(t2)+1], [36.5 36.5],'color',[0 0 0],'linestyle','--');
ylabel('# Spokes Detected'); title('Low Contrast Object Detectability');
xlim([min(t2)-1 max(t2)+1]); 
ylim([min(0.98*min(LCOD),0.90*37) 40])
set(gca,'FontSize',fs); grid off;

subplot(3,3,9);
[t1_sorted,t1_order]=sort(t1);
SNR_sorted = SNR(t1_order);
plot(t1_sorted,SNR_sorted,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% for i=1:length(Dist)
%     plot(t1_sorted(i),SNR_sorted(i),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'MarkerSize',ms); hold on;
% end
s = movingstd(SNR_sorted,k,[]);
M = conv(SNR_sorted, ones(k,1)/k, 'same');
s(1:7)=s(8);s(end-7:end)=s(end-7);
M(1:7)=M(8);M(end-7:end)=M(end-7);
plot(t1_sorted,M,'k-');
plot(t1_sorted,M-2*s,'k--');
plot(t1_sorted,M+2*s,'k--');
ylabel('SNR'); title('Signal To Noise Ratio');
xlim([min(t2)-1 max(t2)+1]); 
%ylim([0 max(1.1*max(Dist),5)]);
set(gca,'FontSize',fs); grid off;

saveas(gcf,[RESULTS_dir filesep 'AllTimeT1.fig'])
print_current_figure(400,[RESULTS_dir filesep 'AllTimeT1.png'])
pause(2); close all;

%% EPI long term

clear T1100Index T1_binary
% grab all with s_all.type = 'T1'
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

figure; set(gcf, 'units','normalized','outerposition',[0 0 1 1]); 
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

%% EPI Single Time Point

name = s_all(EPIMostRecentIndex).name;
clear s
load([OUT_dir filesep name filesep 'Results.mat']) % load S

D = imread([OUT_dir filesep name filesep 'Ghosting.png']);
E = imread([OUT_dir filesep name filesep 'GeometricDistortion.png']);

clear data;
data = {'Scan Name',s.name,[],''; ...
    'Scan Type',s.type ,[],''; ...
    'Scanner',s.scanner,[], ''; ...
    'Date',s.scan_date, [],''; ...
    'Description',s.series_description,[], ''; ...
    'Resolution (xy)',num2str(s.reso(1)),[], '';}   ;   

% Shewhart Chart
% get last 20 (21:end-2) measures (highest t)
[t1_sorted,t1_order]=sort(t1);

Param = GA; rownum = 7; name = 'Geometric Accuracy'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = Dist; rownum = 8; name = 'Distortion'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = PIU; rownum = 9; name = 'Image Intensity Uniformity'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = PSG; rownum = 10; name = 'Ghosting Ratio'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = SNR; rownum = 11; name = 'Signal To Noise Ratio'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = SFNR; rownum = 12; name = 'Signal To Fluctuation Noise Ratio'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = PF; rownum = 13; name = 'Percent Fluctuation'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = Drift; rownum = 14; name = 'Drift'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

Param = RDC; rownum = 15; name = 'Radius of Decorrelation'
Val_sorted = Param(t1_order);
[Shewart,message] = ShewartEval(Val_sorted);
data{rownum,1} = name; data{rownum,2}=Param(end); data{rownum,3}=Shewart; data{rownum,4}=message;

close all;
f = figure(5); set(gcf, 'units','normalized','outerposition',[0 0 1 1]); 
t = uitable(f,'Data',data,'ColumnWidth',{170,110,50,280})
t.ColumnName = {'Measurement','Value','Shewart','Actionable?'};
subplot(1,2,2),plot(3)
pos = get(subplot(1,2,2),'position');
delete(subplot(2,1,2))
set(t,'units','normalized')
set(t,'position',pos)
axis off; 

subplot(2,2,1); imagesc(D); axis equal; axis off; title('Ghosting ROIs');
subplot(2,2,3); imagesc(E); axis equal; axis off; title('Distortion Measures');

saveas(gcf,[RESULTS_dir filesep 'SingleTimeEPI.fig'])
print_current_figure(400,[RESULTS_dir filesep 'SingleTimeEPI.png'])
pause(2); close all;

%% EPI All Time

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

figure; set(gcf, 'units','normalized','outerposition',[0 0 1 1]); 
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
saveas(gcf,[RESULTS_dir filesep 'AllTimeEPI.fig'])
set(gcf, 'PaperPositionMode', 'auto');
print_current_figure(400,[RESULTS_dir filesep 'AllTimeEPI.png'])
pause(2); close all;




