function [distance_real_hori,distance_real_vert] = GeometricAccuracy_S1(img,reso)

I = squeeze(img(:,:,1));
pxl_sz = reso(:);
per = .1;
% Find the water intensity peak on the image intensity histogram. 
% This intensity peak will be used to threshold image.

%1.specify intensity range
I_max=double(max(max(I)));
[hist_cnt,hist_int]=hist(I(:),0:I_max);
hist_sample_start=round(I_max*per);%percentage of max as min to exclude air intenisty
%2.find max within specified histogram range
[int_cnt,int_pk]=max(hist_cnt(hist_sample_start:size(hist_cnt,2)));%v2
int_pk=int_pk+hist_sample_start-1;%v2

mu = int_pk;
% threshold
I_bin = I > mu/2;

band_per=[.3 .6];
% sum each colum position across central rows, get low (L) and high (R) indices
row_band=zeros(1,size(I_bin,2));
for i=round(band_per(1,1)*size(I_bin,1)): round(band_per(1,2)*size(I_bin,1))
    row_band=row_band+double(I_bin(i,:));
end
ind_l=find(row_band>0,1);%v3
ind_r=find(row_band>0,1,'last');%v3

% do same across columns
col_band=zeros(size(I_bin,1),1);
for i=round(band_per(1,1)*size(I_bin,2)):...
        round(band_per(1,2)*size(I_bin,2))
    col_band=col_band+double(I_bin(:,i));
end
ind_t=find(col_band>0,1);%v3
ind_b=find(col_band>0,1,'last');%v3

% find pixel distances and convert to mm
distance_row=ind_r-ind_l;
distance_col=ind_b-ind_t;

distance_real_r=distance_row*pxl_sz(2,1);
distance_real_hori=round(distance_real_r*10)/10;
distance_real_c=distance_col*pxl_sz(2,1);
distance_real_vert=round(distance_real_c*10)/10;

% plot
ind_l_y=(ind_b-ind_t)/2+ind_t;%v5
ind_t_x=(ind_r-ind_l)/2+ind_l;%v5
figure; imshow(I,[]); hold on
plot([ind_t_x,ind_t_x],[ind_t,ind_b],'Color','r','LineWidth',2);%v5
plot([ind_l,ind_r],[ind_l_y,ind_l_y],'Color','r','LineWidth',2);%v5
text(ind_t_x-50,ind_t+20,...%HW:50 pxls to left
    [num2str(distance_real_vert) '\rightarrow'],...%use 'normalized' to scale
    'Color','r','FontUnits','normalized');      %letter to image
text(ind_t_x+20,ind_l_y-20,...%HW:20 pxls up
    ['\downarrow' num2str(distance_real_hori)],...%use 'normalized' to scale
    'Color','r','FontUnits','normalized');     %letter to image
hold off
end

