function [truth,spokes] = LowContrastDetectability(img,reso,slice)
 
    
 % get slice 8-11
Slice = squeeze(img(:,:,slice));
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

  %5.region growth start from phantom centre to make disk mask
  
tolerance=1;
Igray=I_bin;
x=ind_centre(1,1);
y=ind_centre(1,2);
visual=0;

Phi = false(size(Igray,1),size(Igray,2));
ref = true(size(Igray,1),size(Igray,2));
PhiOld = Phi;
Phi(uint8(x),uint8(y)) = 1;
while(sum(Phi(:)) ~= sum(PhiOld(:)))
    PhiOld = Phi;
    segm_val = Igray(Phi);
    meanSeg = mean(segm_val);
    posVoisinsPhi = imdilate(Phi,strel('disk',1,0)) - Phi;
    voisins = find(posVoisinsPhi);
    valeursVoisins = Igray(voisins);
    Phi(voisins(valeursVoisins > ...
        meanSeg - tolerance & valeursVoisins < ...
        meanSeg + tolerance)) = 1;
end

I_mask = Phi;

% find center of disk
band_per=[.3 .6];
% sum each colum position across central rows, get low (L) and high (R) indices
row_band=zeros(1,size(I_mask,2));
for i=round(band_per(1,1)*size(I_mask,1)): round(band_per(1,2)*size(I_mask,1))
    row_band=row_band+double(I_mask(i,:));
end
ind_row_low_d=find(row_band>0,1);%v3
ind_row_high_d=find(row_band>0,1,'last');%v3

% do same across columns
col_band=zeros(size(I_mask,1),1);
for i=round(band_per(1,1)*size(I_mask,2)):...
        round(band_per(1,2)*size(I_mask,2))
    col_band=col_band+double(I_mask(:,i));
end
ind_col_low_d=find(col_band>0,1);%v3
ind_col_high_d=find(col_band>0,1,'last');%v3

ind_centre_d=[round((ind_col_high_d-ind_col_low_d)/2+ind_col_low_d) round((ind_row_high_d-ind_row_low_d)/2+ind_row_low_d)];%centre of disk
      
%7.create line segment, radiating from disk centre
I_masked = I .* I_mask;



radius=round((ind_row_high_d-ind_row_low_d)/2);%radius of disk

%8.find radius of 3 spokes
r_mm_spk1=12.6;%HW:mean of measurement on S11
r_mm_spk2=25.3;%HW:mean of measurement on S11
r_mm_spk3=38.1;%HW:mean of measurement on S11

r_pxl_spk1=round(r_mm_spk1/pxl_sz(1,1));
r_pxl_spk2=round(r_mm_spk2/pxl_sz(1,1));
r_pxl_spk3=round(r_mm_spk3/pxl_sz(1,1));

 %9.sample circumference intensity of 3 radius for spokes
 r_pxl=r_pxl_spk1;
 rad_interval=.01;
 % define the centre of circle & angle interval and find intensity along circumference
x0=ind_centre_d(1,2);
y0=ind_centre_d(1,1);
theta=0:rad_interval:2*pi;
Int_circum=zeros(1,1);
cnt=1;
for i=1:size(theta,2)
    xi=r_pxl*cos(3*pi/2+theta(1,i))+x0; yi=r_pxl*sin(3*pi/2+theta(1,i))+y0;
    Int_circum(1,cnt)=I_masked(round(yi),round(xi));
    cnt=cnt+1;
end
 Int_circum_1 = Int_circum;
 
 r_pxl=r_pxl_spk2;
 rad_interval=.01;
Int_circum=zeros(1,1);
cnt=1;
for i=1:size(theta,2)
    xi=r_pxl*cos(3*pi/2+theta(1,i))+x0; yi=r_pxl*sin(3*pi/2+theta(1,i))+y0;
    Int_circum(1,cnt)=I_masked(round(yi),round(xi));
    cnt=cnt+1;
end
 Int_circum_2 = Int_circum;
 
  r_pxl=r_pxl_spk3;
  rad_interval=.01;
Int_circum=zeros(1,1);
cnt=1;
for i=1:size(theta,2)
    xi=r_pxl*cos(3*pi/2+theta(1,i))+x0; yi=r_pxl*sin(3*pi/2+theta(1,i))+y0;
    Int_circum(1,cnt)=I_masked(round(yi),round(xi));
    cnt=cnt+1;
end
 Int_circum_3 = Int_circum;
 
 found = zeros(10,3);
 
 min_distance = 24; % degrees
 min_distance_points = min_distance*length(Int_circum_1)/360;
 Biggest_kernel1 = 1.5*7/(2*pi*r_mm_spk1) * length(Int_circum_1);
 Smallest_kernel1 = 1.5*1.5/(2*pi*r_mm_spk1) * length(Int_circum_1);
 mid1 = Biggest_kernel1 - (Biggest_kernel1-Smallest_kernel1)/3;
 mid2 = Biggest_kernel1 - 2*(Biggest_kernel1-Smallest_kernel1)/3;
 Int1 = [Int_circum_1, Int_circum_1, Int_circum_1];
 Int1_smooth = smooth(Int1,Biggest_kernel1,'sgolay') + smooth(Int1,mid1,'sgolay') + smooth(Int1,mid2,'sgolay') + smooth(Int1,Smallest_kernel1,'sgolay');
 One =  smooth(Int1_smooth,min_distance_points/2);
 [pks,locs] = findpeaks(One,1:length(One),'MinPeakDistance',min_distance_points);
locs(locs<length(Int_circum_1)+1 | locs>length(Int_circum_1)+length(Int_circum_1)) = [];
locs = locs - length(Int_circum_1);
locs_bin = ceil(locs/(length(Int_circum_1)/10));

for i=1:length(locs_bin)
    
    found(locs_bin(i),1) = 1;
end

%figure; plot(Int1); figure; plot(Int1_smooth); figure; plot(One);
 min_distance_points = min_distance*length(Int_circum_2)/360;
 Biggest_kernel1 = 1.5*7/(2*pi*r_mm_spk2) * length(Int_circum_2);
 Smallest_kernel1 = 1.5*1.5/(2*pi*r_mm_spk2) * length(Int_circum_2);
 mid1 = Biggest_kernel1 - (Biggest_kernel1-Smallest_kernel1)/3;
 mid2 = Biggest_kernel1 - 2*(Biggest_kernel1-Smallest_kernel1)/3;
 Int1 = [Int_circum_2, Int_circum_2, Int_circum_2];
 Int1_smooth = smooth(Int1,Biggest_kernel1,'sgolay') + smooth(Int1,mid1,'sgolay') + smooth(Int1,mid2,'sgolay') + smooth(Int1,Smallest_kernel1,'sgolay');
 One =  smooth(Int1_smooth,min_distance_points/2);
 [pks,locs] = findpeaks(One,1:length(One),'MinPeakDistance',min_distance_points);
locs(locs<length(Int_circum_2)+1 | locs>length(Int_circum_2)+length(Int_circum_2)) = [];
locs = locs - length(Int_circum_2);
locs_bin = ceil(locs/(length(Int_circum_2)/10));

for i=1:length(locs_bin)
    
    found(locs_bin(i),2) = 1;
end

  min_distance_points = min_distance*length(Int_circum_3)/360;
 Biggest_kernel1 = 1.5*7/(2*pi*r_mm_spk3) * length(Int_circum_3);
 Smallest_kernel1 = 1.5*1.5/(2*pi*r_mm_spk3) * length(Int_circum_3);
 mid1 = Biggest_kernel1 - (Biggest_kernel1-Smallest_kernel1)/3;
 mid2 = Biggest_kernel1 - 2*(Biggest_kernel1-Smallest_kernel1)/3;
 Int1 = [Int_circum_2, Int_circum_2, Int_circum_2];
 Int1_smooth = smooth(Int1,Biggest_kernel1,'sgolay') + smooth(Int1,mid1,'sgolay') + smooth(Int1,mid2,'sgolay') + smooth(Int1,Smallest_kernel1,'sgolay');
 One =  smooth(Int1_smooth,min_distance_points/2);
 [pks,locs] = findpeaks(One,1:length(One),'MinPeakDistance',min_distance_points);
locs(locs<length(Int_circum_3)+1 | locs>length(Int_circum_3)+length(Int_circum_3)) = [];
locs = locs - length(Int_circum_3);
locs_bin = ceil(locs/(length(Int_circum_3)/10));
 
 
for i=1:length(locs_bin)
    
    found(locs_bin(i),3) = 1;
end

 
 truth = min(found,[],2);
 spokes = sum(truth);
 
 figure; imagesc(I_masked); colormap gray;
 
 
 
 
%  %10.remove repeated intensities from profile
%  cnt=1;
%  Int_circum_1_f=0;
%  for j=1:size(Int_circum_1,2)-1
%      if Int_circum_1(1,j+1)~=Int_circum_1(1,j)
%          Int_circum_1_f(1,cnt)=Int_circum_1(1,j);
%          cnt=cnt+1;
%      end
%  end
%  cnt=1;
%  Int_circum_2_f=0;
%  for j=1:size(Int_circum_2,2)-1
%      if Int_circum_2(1,j+1)~=Int_circum_2(1,j)
%          Int_circum_2_f(1,cnt)=Int_circum_2(1,j);
%          cnt=cnt+1;
%      end
%  end
%  cnt=1;
%  Int_circum_3_f=0;
%  for j=1:size(Int_circum_3,2)-1
%      if Int_circum_3(1,j+1)~=Int_circum_3(1,j)
%          Int_circum_3_f(1,cnt)=Int_circum_3(1,j);
%          cnt=cnt+1;
%      end
%  end
 
%  %10.calculate radius to the inner edge of 3 spokes with largest radius
%  spoke_dia_mm=7;
%  r_pxl_inner1=r_pxl_spk1-round(spoke_dia_mm/(2*pxl_sz(1,1)));
%  r_pxl_inner2=r_pxl_spk2-round(spoke_dia_mm/(2*pxl_sz(1,1)));
%  r_pxl_inner3=r_pxl_spk3-round(spoke_dia_mm/(2*pxl_sz(1,1)));
%  
%  %11.sample circumference intensity of 3 radius for background
%  r_pxl=r_pxl_inner1;
%  rad_interval=0.16;
%  x0=ind_centre_d(1,2);
%  y0=ind_centre_d(1,1);
%  theta=0:rad_interval:2*pi;
%  Int_circum=zeros(1,1);
%  cnt=1;
%  for i=1:size(theta,2)
%      xi=r_pxl*cos(3*pi/2+theta(1,i))+x0;
%      yi=r_pxl*sin(3*pi/2+theta(1,i))+y0;
%      Int_circum(1,cnt)=I_masked(round(yi),round(xi));
%      cnt=cnt+1;
%  end
%  Int_circum_bkgd_1 = Int_circum;
%  
%  r_pxl=r_pxl_inner2;
%  rad_interval=0.1;
%  Int_circum=zeros(1,1);
%  cnt=1;
%  for i=1:size(theta,2)
%      xi=r_pxl*cos(3*pi/2+theta(1,i))+x0;
%      yi=r_pxl*sin(3*pi/2+theta(1,i))+y0;
%      Int_circum(1,cnt)=I_masked(round(yi),round(xi));
%      cnt=cnt+1;
%  end
%  Int_circum_bkgd_2 = Int_circum;
%  
%  r_pxl=r_pxl_inner3;
%  rad_interval=0.1;
%  Int_circum=zeros(1,1);
%  cnt=1;
%  for i=1:size(theta,2)
%      xi=r_pxl*cos(3*pi/2+theta(1,i))+x0;
%      yi=r_pxl*sin(3*pi/2+theta(1,i))+y0;
%      Int_circum(1,cnt)=I_masked(round(yi),round(xi));
%      cnt=cnt+1;
%  end
%  Int_circum_bkgd_3 = Int_circum;
%  
%  r_pxl=r_pxl_inner1-1;
%  rad_interval=0.16;
%  Int_circum=zeros(1,1);
%  cnt=1;
%  for i=1:size(theta,2)
%      xi=r_pxl*cos(3*pi/2+theta(1,i))+x0;
%      yi=r_pxl*sin(3*pi/2+theta(1,i))+y0;
%      Int_circum(1,cnt)=I_masked(round(yi),round(xi));
%      cnt=cnt+1;
%  end
%  Int_circum_bkgd_1_dash = Int_circum;
%  
%  r_pxl=r_pxl_inner2-1;
%  rad_interval=0.1;
%  Int_circum=zeros(1,1);
%  cnt=1;
%  for i=1:size(theta,2)
%      xi=r_pxl*cos(3*pi/2+theta(1,i))+x0;
%      yi=r_pxl*sin(3*pi/2+theta(1,i))+y0;
%      Int_circum(1,cnt)=I_masked(round(yi),round(xi));
%      cnt=cnt+1;
%  end
%  Int_circum_bkgd_2_dash = Int_circum;
%  
% r_pxl=r_pxl_inner3-1;
%  rad_interval=0.1;
%  Int_circum=zeros(1,1);
%  cnt=1;
%  for i=1:size(theta,2)
%      xi=r_pxl*cos(3*pi/2+theta(1,i))+x0;
%      yi=r_pxl*sin(3*pi/2+theta(1,i))+y0;
%      Int_circum(1,cnt)=I_masked(round(yi),round(xi));
%      cnt=cnt+1;
%  end
%  Int_circum_bkgd_3_dash = Int_circum;
%  
%  %12.find the normalised spoke intensity
%  spk1_norm=Int_circum_1./Int_circum_bkgd_1;
%  spk2_norm=Int_circum_2./Int_circum_bkgd_2;
%  spk3_norm=Int_circum_3./Int_circum_bkgd_3;
%  
%  %13.find the normalised background intensity
%  bkgd1_norm=Int_circum_bkgd_1./Int_circum_bkgd_1_dash;
%  bkgd2_norm=Int_circum_bkgd_2./Int_circum_bkgd_2_dash;
%  bkgd3_norm=Int_circum_bkgd_3./Int_circum_bkgd_3_dash;
%  
%          %14.find the spoke peaks
%         [spk1_pks,spk1_pks_loc]=fun_ACR_FindSpokePeak...
%             ([],1,spk1_norm,0.16,0);
%         [spk2_pks,spk2_pks_loc]=fun_ACR_FindSpokePeak...
%             ([],2,spk2_norm,0.1,0);
%         [spk3_pks,spk3_pks_loc]=fun_ACR_FindSpokePeak...
%             ([],3,spk3_norm,0.1,0);
% 


%% 


     



%         %14.find the spoke peaks
%         [spk1_pks,spk1_pks_loc]=fun_ACR_FindSpokePeak...
%             ([],1,spk1_norm,0.16,visual);
%         [spk2_pks,spk2_pks_loc]=fun_ACR_FindSpokePeak...
%             ([],2,spk2_norm,0.1,visual);
%         [spk3_pks,spk3_pks_loc]=fun_ACR_FindSpokePeak...
%             ([],3,spk3_norm,0.1,visual);
%         %15.find the min values between 2 peaks
%         [spk1_vly,spk1_vly_loc]=fun_ACR_FindSpokeValley...
%             ([],1,spk1_norm,spk1_pks_loc,visual);
%         [spk2_vly,spk2_vly_loc]=fun_ACR_FindSpokeValley...
%             ([],2,spk2_norm,spk2_pks_loc,visual);
%         [spk3_vly,spk3_vly_loc]=fun_ACR_FindSpokeValley...
%             ([],3,spk3_norm,spk3_pks_loc,visual);
%         %16.use peak and min to find the contrast
%         spk1_vly_diff=spk1_pks-spk1_vly;
%         spk2_vly_diff=spk2_pks-spk2_vly;
%         spk3_vly_diff=spk3_pks-spk3_vly;
%         %17.compare the calculated contrast to given contrast
%         dummy_diff=cat(1,spk1_vly_diff,spk2_vly_diff,spk3_vly_diff);
%         dummy_bkgrd=cat(1,std(bkgd1_norm),std(bkgd2_norm),std(bkgd3_norm));
%         disp('LCOD shows if can identify spoke on r=1,2,3:');
%         for i=1:10
%             LCOD(:,i)=dummy_diff(:,i)>dummy_bkgrd;
%         end
%         %18.count the spokes
%         cnt=0;
%         for i=1:10
%             if sum(LCOD(:,i))==3
%                 cnt=cnt+1;
%             else
%                 break;
%             end
%         end
        %bug:if cannot find spoke peak, function stops (S8&S9)
        %1.try -2 radius for background or half radius for background
        %2.add +/- 1 pxl radius uncertainty to all 3 radius of spoke intensity
        %3.can i use the normalised infor or other info to find the contrast value
        %of the current disk (yes i can and the value is < expected value)
      
end
