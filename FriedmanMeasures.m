function [SFNR,SNR,PF,Drift,RDC]=FriedmanMeasures(img,reso,slice)
%ROI: 1=water ROI, 2=top ROI, 3=bottom ROI, 4=left ROI,5=right ROI

Slice = squeeze(img(:,:,slice,round(size(img,4)/2)));

% get slice 7
%Slice = squeeze(img(:,:,7));
pxl_sz = reso(:);
I = imrotate(flip(Slice,1),90);
per = .1;
I_max=double(max(max(I)));
[hist_cnt,hist_int]=hist(I(:),0:I_max);
hist_sample_start=round(I_max*per);%percentage of max as min to exclude air intenisty
%2.find max within specified histogram range
[int_cnt,int_pk]=max(hist_cnt(hist_sample_start:size(hist_cnt,2)));%v2
int_pk=int_pk+hist_sample_start-1;%v2

water_mean=int_pk;
I_bin = I > water_mean/3;

% find center of phantom
band_per=[.3 .6];
% sum each colum position across central rows, get low (L) and high (R) indices
row_band=zeros(1,size(I_bin,2));
for i=round(band_per(1,1)*size(I_bin,1)): round(band_per(1,2)*size(I_bin,1))
    row_band=row_band+double(I_bin(i,:));
end
ind_row_low=find(row_band>0,1);%v3
ind_row_high=find(row_band>0,1,'last');%v3

% do same across columns
col_band=zeros(size(I_bin,1),1);
for i=round(band_per(1,1)*size(I_bin,2)):...
        round(band_per(1,2)*size(I_bin,2))
    col_band=col_band+double(I_bin(:,i));
end
ind_col_low=find(col_band>0,1);%v3
ind_col_high=find(col_band>0,1,'last');%v3

ind_centre=[round((ind_col_high-ind_col_low)/2+ind_col_low) round((ind_row_high-ind_row_low)/2+ind_row_low)];%phantom centre (y,x)

%6.create circular ROI around phantom centre with user specified radius
 area_ROI_water=190;
 
 radius_cm=sqrt(area_ROI_water/pi);
radius_mm=radius_cm*10;
radius_pxl=radius_mm/pxl_sz(1,1);%from this line, use ellipse ROI if image
theta=0:0.01:2*pi;               %is anisotropic (see test.m for working)
x=radius_pxl*cos(theta)+ind_centre(1,2);
y=radius_pxl*sin(theta)+ind_centre(1,1);%HW:centre of ROI is 10 pxls%below phantom centre

BW=roipoly(I,x,y);

I_ROI_water = I.*BW;
I_ROI_water_sum=sum(I(BW));%this finds sum of water within mask
I_ROI_water_mean=I_ROI_water_sum/size(I(BW),1);%mean of water ROI
I_ROI_water_sigma=std(double(I(BW)));%std of water ROI
mu_S7=I_ROI_water_mean;

%% SIGNAL IMAGE
% simple average, voxel by voxel, across 198 images
img_rotate = imrotate(flip(img,1),90);
img_end = squeeze(img_rotate(:,:,slice,3:end));
SignalImage = mean(img_end,3);

%% TEMPORAL FLUCTUATION NOISE IMAGE
% time series across 198 images for each voxel is detremended with 2nd
% order polynomial
% tfni is SD of residuals, voxel by voxel, after detrending

tfni = zeros(size(img_end,1),size(img_end,2));
t = 1:size(img_end,3); t = t(:);
for rows = 1:size(img_end,1)
    for cols = 1:size(img_end,2)
        if I_bin(rows,cols)==1
            y = squeeze(img_end(rows,cols,:));
            p = polyfit(t,y,2);
            y2 = polyval(p,t);
            
            %figure
            %plot(t,y,'o',t,y2)
            %title('Plot of Data (Points) and Model (Line)')
            res = y - y2;
            
            %figure, plot(t,res,'+')
            %title('Plot of the Residuals')
            
            tfni(rows,cols) = std(res);
        end
    end
end

%% SIGNAL TO FLUCTUATION NOISE RATIO IMAGE
% signal image and tfi divded voxel by voxel to create SFNR, a 21by21 voxel
% ROI, pl;aced in center of image, is created, average of SFNR across 441
% voxels is SFNR summary
        
SFNRI = SignalImage./tfni;

% 21 by 21 region centered on image
% ind_centre is row, col
ROI = zeros(size(SFNRI));
ROI(ind_centre(1)-10:ind_centre(1)+10,ind_centre(2)-10:ind_centre(2)+10) = 1;
SFNR = mean(SFNRI(ROI==1));

%% STATIC SPATIAL NOISE IMAGE
% sum all odd number images
% sum all even numbered images
% difference betw3en sums is taken as raw measure of static spatial noise

numImages = size(img_end,3);
sumODD = sum(img_end(:,:,1:2:end),3);
sumEVEN = sum(img_end(:,:,2:2:end),3);
DIFF = sumODD-sumEVEN;

%% SNR Summary value
% static spatial noise variance summary value is variance of static spatial
% noise (DIFF) image across 21by21 voxel ROI centered on image
% signal summary value is average of signal image across same ROI
% SNR = (signalsummaryvalue)/sqrt((variancesummaryvalue)/XtimePoints)

ssnvs = var(DIFF(ROI==1));
signalsummaryvalue = mean(SignalImage(ROI==1));
SNR = (signalsummaryvalue)/sqrt((ssnvs)/numImages);

%% PERCENT FLUCTUATION
% time series of average intensity of 21by21 ROI centered on image is
% obtained
% second order polynomial trend fit to data
% mean signal intensity of time series (prior to detrend) and SD of
% residuals after subtracting fit line from data computed
% PF = 100*(SDofResiduals)/(meanSignalIntensity)

timeseries=zeros(numImages,1);
for mm = 1:numImages
    curr_image = squeeze(img_end(:,:,mm));
    curr_signal = mean(curr_image(ROI==1));
    timeseries(mm) = curr_signal;
end

t = 1:size(img_end,3); t = t(:);

p = polyfit(t,timeseries,2);
y2 = polyval(p,t);
%figure; plot(t,timeseries,'o',t,y2); title('Plot of Data (Points) and Model (Line)')
resids = timeseries - y2;
%figure, plot(t,res,'+'); title('Plot of the Residuals')
meanSignalIntensityTimeSeries = mean(timeseries);
SDresid = std(resids);
PF = 100*(SDresid)/(meanSignalIntensityTimeSeries);

%% DRIFT
% subtracting the minimum fit value from the maximum fit and dividing by
% mean signal intensity
% drift value multiplied by 100 to obtain percentage

Drift = 100*(max(y2)-min(y2))/meanSignalIntensityTimeSeries;

%% 8)Weisskoff Analysis (Friedman, JMRI)/ Radius of Decorrelation


% calculate CV of time series (SD of time series divided by mean of time
% series) for ROIs of differing sizes (from 1 to 21, odd)
ind=1;
CV = zeros(length(1:2:21),1);
t = 1:size(img_end,3); t = t(:);
for ROIsize = 1:2:21
    % get ROI
    ROI_curr = zeros(size(SFNRI));
    ROI_curr(ind_centre(1)-floor(ROIsize/2):ind_centre(1)+floor(ROIsize/2),ind_centre(2)-floor(ROIsize/2):ind_centre(2)+floor(ROIsize/2)) = 1;
    % get time series
    for mm = 1:numImages
        curr_image = squeeze(img_end(:,:,mm));
        curr_signal = mean(curr_image(ROI_curr==1));
        timeseries(mm) = curr_signal;
    end
    % detrend time series
    

p = polyfit(t,timeseries,2);
y2 = polyval(p,t);
%figure; plot(t,timeseries,'o',t,y2); title('Plot of Data (Points) and Model (Line)')
resids = timeseries - y2;
    
    CV(ind) = 100*std(resids)/mean(timeseries);
    ind=ind+1;
end

 RDC = CV(1)/CV(end);
 % RDC may be thought of as measures of size of ROI at which statistical
 % independence of voxels is lost

end










