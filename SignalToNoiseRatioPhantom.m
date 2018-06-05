function  [SNR] = SignalToNoiseRatioPhantom(img,reso)

% get slice 7, get ROI mostly covering circle
% get background of something in top right corner!
Slice = squeeze(img(:,:,7));
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
I_bin = I > water_mean/2;

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
 area_ROI_water=200;
 
 radius_cm=sqrt(area_ROI_water/pi);
radius_mm=radius_cm*10;
radius_pxl=radius_mm/pxl_sz(1,1);%from this line, use ellipse ROI if image
theta=0:0.01:2*pi;               %is anisotropic (see test.m for working)
x=radius_pxl*cos(theta)+ind_centre(1,2);
y=radius_pxl*sin(theta)+ind_centre(1,1)+10;%HW:centre of ROI is 10 pxls%below phantom centre

BW=roipoly(I,x,y);

I_ROI_water = I.*BW;
I_ROI_water_sum=sum(I(BW));%this finds sum of water within mask
I_ROI_water_mean=I_ROI_water_sum/size(I(BW),1);%mean of water ROI
I_ROI_water_sigma=std(double(I(BW)));%std of water ROI
mu_S7=I_ROI_water_mean;


% get top index (lowest row) and right index (highest column) of I_bin
sumAcrossRows = sum(I_bin,2);
RowLowIndex=find(sumAcrossRows>0,1)-1; % last row that has all zeros
sumAcrossCols = sum(I_bin,1);
ColHighIndex=find(col_band>0,1,'last')+1; % first column that has all zeros

% find number of cols, subtract 10% from each side, get indices of rows
rowIndices = [1+round(length(1:RowLowIndex)*.1) RowLowIndex-round(length(1:RowLowIndex)*.1)];

colIndices =  [ColHighIndex + round(length(ColHighIndex:size(I_bin,2))*.1) size(I_bin,2)-round(length(ColHighIndex:size(I_bin,2))*.1)];

bckground = I(rowIndices(1):rowIndices(2),colIndices(1):colIndices(2));
stdbg = std(bckground(:));

SNR = sqrt(2*(1-pi/4))    *   mu_S7/stdbg;

% figure of image with ROI and mask outlined
bgm = zeros(size(BW));
bgm(rowIndices(1):rowIndices(2),colIndices(1):colIndices(2))=1;

B = bwboundaries(BW);
C = bwboundaries(bgm);

figure; imagesc(I); axis equal; colormap gray; axis off;
hold on;
b = B{1};
c = C{1};
plot(b(:,2),b(:,1),'g','linewidth',2);
plot(c(:,2),c(:,1),'r','linewidth',2);



end