function [] = T1_Axial_Eval(s_all,Dates,RESULTS_dir)

T1_binary = strcmp({s_all.type}, 'T1'); % grab all with s_all.type = 'T1'
T1Index = find(T1_binary);
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

figure; set(gcf, 'units','normalized','outerposition',[0.05 0.05 .9 .9]); 
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
textDist = t1(Dist_sorted<M-2*s | Dist_sorted>M+2*s);

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
textHCSR = t1(HCSR>1.05);

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
textST = t1(SliceThick<4.3 | SliceThick>5.7);

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
textSP = t1(SlicePos<-2.5 | SlicePos>2.5);

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
textPIU = t1(PIU<.82);

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
textG = t1(PSG>.025);

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
textLC = t1(LCOD<36.5);

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
textSNR = t1(SNR<M-2*s | SNR>M+2*s);

%% Save Figures
% All Time, 100 days, 30 days, 7 days

drawnow; 
pause(1);
saveas(gcf,[RESULTS_dir filesep 'T1_AllTime.fig'])
pause(1);
print_current_figure(400,[RESULTS_dir filesep 'T1_AllTime.png'])
pause(1);

% All Time
subplot(3,3,1); xlim([max(t2)-99 max(t2)]); 
subplot(3,3,2); xlim([max(t2)-99 max(t2)]); 
subplot(3,3,3); xlim([max(t2)-99 max(t2)]); 
subplot(3,3,4); xlim([max(t2)-99 max(t2)]); 
subplot(3,3,5); xlim([max(t2)-99 max(t2)]); 
subplot(3,3,6); xlim([max(t2)-99 max(t2)]); 
subplot(3,3,7); xlim([max(t2)-99 max(t2)]); 
subplot(3,3,8); xlim([max(t2)-99 max(t2)]); 
subplot(3,3,9); xlim([max(t2)-99 max(t2)]); 

drawnow; 
pause(1);
saveas(gcf,[RESULTS_dir filesep 'T1_100Days.fig'])
pause(1);
print_current_figure(400,[RESULTS_dir filesep 'T1_100Days.png'])
pause(1);

subplot(3,3,1); xlim([max(t2)-30 max(t2)]); 
subplot(3,3,2); xlim([max(t2)-30 max(t2)]); 
subplot(3,3,3); xlim([max(t2)-30 max(t2)]); 
subplot(3,3,4); xlim([max(t2)-30 max(t2)]); 
subplot(3,3,5); xlim([max(t2)-30 max(t2)]); 
subplot(3,3,6); xlim([max(t2)-30 max(t2)]); 
subplot(3,3,7); xlim([max(t2)-30 max(t2)]); 
subplot(3,3,8); xlim([max(t2)-30 max(t2)]); 
subplot(3,3,9); xlim([max(t2)-30 max(t2)]); 

drawnow; 
pause(1);
saveas(gcf,[RESULTS_dir filesep 'T1_30Days.fig'])
pause(1);
print_current_figure(400,[RESULTS_dir filesep 'T1_30Days.png'])
pause(1);

subplot(3,3,1); xlim([max(t2)-6 max(t2)]); 
subplot(3,3,2); xlim([max(t2)-6 max(t2)]); 
subplot(3,3,3); xlim([max(t2)-6 max(t2)]); 
subplot(3,3,4); xlim([max(t2)-6 max(t2)]); 
subplot(3,3,5); xlim([max(t2)-6 max(t2)]); 
subplot(3,3,6); xlim([max(t2)-6 max(t2)]); 
subplot(3,3,7); xlim([max(t2)-6 max(t2)]); 
subplot(3,3,8); xlim([max(t2)-6 max(t2)]); 
subplot(3,3,9); xlim([max(t2)-6 max(t2)]); 

drawnow; 
pause(1);
saveas(gcf,[RESULTS_dir filesep 'T1_7Days.fig'])
pause(1);
print_current_figure(400,[RESULTS_dir filesep 'T1_7Days.png'])
pause(1);
close all

%% Text - All Time

fileID = fopen([RESULTS_dir filesep 'T1_AllTime.txt'],'w');

fprintf(fileID,'Geometric Accuracy Failures:\n');
for i =1:length(textGA)
    fprintf(fileID,[datestr(textGA(i)) '\n']);
end
fprintf(fileID,'\n');

fprintf(fileID,'Geometric Distortion out of 2 St.Dev. Range:\n');
for i =1:length(textDist)
    fprintf(fileID,[datestr(textDist(i)) '\n']);
end
fprintf(fileID,'\n');

fprintf(fileID,'High Contrast Spatial Res Failures:\n');
for i =1:length(textHCSR)
    fprintf(fileID,[datestr(textHCSR(i)) '\n']);
end
fprintf(fileID,'\n');

fprintf(fileID,'Slice Thickness Failures:\n');
for i =1:length(textST)
    fprintf(fileID,[datestr(textST(i)) '\n']);
end
fprintf(fileID,'\n');

fprintf(fileID,'Slice Position Failures:\n');
for i =1:length(textSP)
    fprintf(fileID,[datestr(textSP(i)) '\n']);
end
fprintf(fileID,'\n');

fprintf(fileID,'Ghosting Ratio Failures:\n');
for i =1:length(textG)
    fprintf(fileID,[datestr(textG(i)) '\n']);
end
fprintf(fileID,'\n');

fprintf(fileID,'SNR out of 2 St.Dev. Range:\n');
for i =1:length(textSNR)
    fprintf(fileID,[datestr(textSNR(i)) '\n']);
end
fprintf(fileID,'\n');

fprintf(fileID,'Low Contrast Detectability Failure:\n');
for i =1:length(textLC)
    fprintf(fileID,[datestr(textLC(i)) '\n']);
end
fprintf(fileID,'\n');

fprintf(fileID,'Image Uniformity Failure:\n');
for i =1:length(textPIU)
    fprintf(fileID,[datestr(textPIU(i)) '\n']);
end
fprintf(fileID,'\n');

fclose(fileID);

%% Text - 100 days

fileID = fopen([RESULTS_dir filesep 'T1_100Days.txt'],'w');

fprintf(fileID,'Geometric Accuracy Failures:\n');
for i =1:length(textGA)
    if datenum(max(t1))-datenum(textGA(i)) < 100
        fprintf(fileID,[datestr(textGA(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Geometric Distortion out of 2 St.Dev. Range:\n');
for i =1:length(textDist)
    if datenum(max(t1))-datenum(textDist(i)) < 100
    fprintf(fileID,[datestr(textDist(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'High Contrast Spatial Res Failures:\n');
for i =1:length(textHCSR)
    if datenum(max(t1))-datenum(textHCSR(i)) < 100
    fprintf(fileID,[datestr(textHCSR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Slice Thickness Failures:\n');
for i =1:length(textST)
    if datenum(max(t1))-datenum(textST(i)) < 100
    fprintf(fileID,[datestr(textST(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Slice Position Failures:\n');
for i =1:length(textSP)
    if datenum(max(t1))-datenum(textSP(i)) < 100
    fprintf(fileID,[datestr(textSP(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Ghosting Ratio Failures:\n');
for i =1:length(textG)
    if datenum(max(t1))-datenum(textG(i)) < 100
    fprintf(fileID,[datestr(textG(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'SNR out of 2 St.Dev. Range:\n');
for i =1:length(textSNR)
    if datenum(max(t1))-datenum(textSNR(i)) < 100
    fprintf(fileID,[datestr(textSNR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Low Contrast Detectability Failure:\n');
for i =1:length(textLC)
    if datenum(max(t1))-datenum(textLC(i)) < 100
    fprintf(fileID,[datestr(textLC(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Image Uniformity Failure:\n');
for i =1:length(textPIU)
    if datenum(max(t1))-datenum(textPIU(i)) < 100
    fprintf(fileID,[datestr(textPIU(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fclose(fileID);

%% Text - 30 days

fileID = fopen([RESULTS_dir filesep 'T1_30Days.txt'],'w');

fprintf(fileID,'Geometric Accuracy Failures:\n');
for i =1:length(textGA)
    if datenum(max(t1))-datenum(textGA(i)) < 30
        fprintf(fileID,[datestr(textGA(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Geometric Distortion out of 2 St.Dev. Range:\n');
for i =1:length(textDist)
    if datenum(max(t1))-datenum(textDist(i)) < 30
    fprintf(fileID,[datestr(textDist(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'High Contrast Spatial Res Failures:\n');
for i =1:length(textHCSR)
    if datenum(max(t1))-datenum(textHCSR(i)) < 30
    fprintf(fileID,[datestr(textHCSR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Slice Thickness Failures:\n');
for i =1:length(textST)
    if datenum(max(t1))-datenum(textST(i)) < 30
    fprintf(fileID,[datestr(textST(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Slice Position Failures:\n');
for i =1:length(textSP)
    if datenum(max(t1))-datenum(textSP(i)) < 30
    fprintf(fileID,[datestr(textSP(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Ghosting Ratio Failures:\n');
for i =1:length(textG)
    if datenum(max(t1))-datenum(textG(i)) < 30
    fprintf(fileID,[datestr(textG(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'SNR out of 2 St.Dev. Range:\n');
for i =1:length(textSNR)
    if datenum(max(t1))-datenum(textSNR(i)) < 30
    fprintf(fileID,[datestr(textSNR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Low Contrast Detectability Failure:\n');
for i =1:length(textLC)
    if datenum(max(t1))-datenum(textLC(i)) < 30
    fprintf(fileID,[datestr(textLC(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Image Uniformity Failure:\n');
for i =1:length(textPIU)
    if datenum(max(t1))-datenum(textPIU(i)) < 30
    fprintf(fileID,[datestr(textPIU(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fclose(fileID);

%% Text - 7 days

fileID = fopen([RESULTS_dir filesep 'T1_7Days.txt'],'w');

fprintf(fileID,'Geometric Accuracy Failures:\n');
for i =1:length(textGA)
    if datenum(max(t1))-datenum(textGA(i)) < 7
        fprintf(fileID,[datestr(textGA(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Geometric Distortion out of 2 St.Dev. Range:\n');
for i =1:length(textDist)
    if datenum(max(t1))-datenum(textDist(i)) < 7
    fprintf(fileID,[datestr(textDist(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'High Contrast Spatial Res Failures:\n');
for i =1:length(textHCSR)
    if datenum(max(t1))-datenum(textHCSR(i)) < 7
    fprintf(fileID,[datestr(textHCSR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Slice Thickness Failures:\n');
for i =1:length(textST)
    if datenum(max(t1))-datenum(textST(i)) < 7
    fprintf(fileID,[datestr(textST(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Slice Position Failures:\n');
for i =1:length(textSP)
    if datenum(max(t1))-datenum(textSP(i)) < 7
    fprintf(fileID,[datestr(textSP(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Ghosting Ratio Failures:\n');
for i =1:length(textG)
    if datenum(max(t1))-datenum(textG(i)) < 7
    fprintf(fileID,[datestr(textG(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'SNR out of 2 St.Dev. Range:\n');
for i =1:length(textSNR)
    if datenum(max(t1))-datenum(textSNR(i)) < 7
    fprintf(fileID,[datestr(textSNR(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Low Contrast Detectability Failure:\n');
for i =1:length(textLC)
    if datenum(max(t1))-datenum(textLC(i)) < 7
    fprintf(fileID,[datestr(textLC(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'Image Uniformity Failure:\n');
for i =1:length(textPIU)
    if datenum(max(t1))-datenum(textPIU(i)) < 7
    fprintf(fileID,[datestr(textPIU(i)) '\n']);
    end
end
fprintf(fileID,'\n');

fclose(fileID);









