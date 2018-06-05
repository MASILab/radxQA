function [distance_real_hori,distance_real_vert,distance_real_ng,distance_real_pg] = GeometricAccuracy_S5(img,reso)

I = squeeze(img(:,:,5));
pxl_sz = reso(:);
per = .1;

%1.specify intensity range
I_max=double(max(max(I)));
[hist_cnt,hist_int]=hist(I(:),0:I_max);
hist_sample_start=round(I_max*per);%percentage of max as min to exclude air intenisty
%2.find max within specified histogram range
[int_cnt,int_pk]=max(hist_cnt(hist_sample_start:size(hist_cnt,2)));%v2
int_pk=int_pk+hist_sample_start-1;%v2

mu = int_pk;

% threshold at mu/2 of half water mean
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

%6.find pixel distance and convert to mm
distance_col=ind_r-ind_l;
distance_row=ind_b-ind_t;
distance_real_c=distance_col*pxl_sz(1,1);
distance_real_hori=round(distance_real_c*10)/10;
distance_real_r=distance_row*pxl_sz(1,1);
distance_real_vert=round(distance_real_r*10)/10;

%7.find rest coord of bndry
ind_l_y=(ind_b-ind_t)/2+ind_t;%v5
ind_t_x=(ind_r-ind_l)/2+ind_l;%v5
figure; subplot(1,2,1); 
imshow(I,[]); hold on

plot([ind_t_x,ind_t_x],[ind_t,ind_b],'Color','r','LineWidth',2);%v5
plot([ind_l,ind_r],[ind_l_y,ind_l_y],'Color','r','LineWidth',2);%v5
text(ind_t_x-50,ind_t+20,...%HW:50 pxls to left
    [num2str(distance_real_vert) '\rightarrow'],...%use 'normalized' to
    'Color','r','FontUnits','normalized');         %scale letter to image
text(ind_t_x+20,ind_l_y-10,...%HW:20 pxls up
    ['\downarrow' num2str(distance_real_hori)],...%use 'normalized' to
    'Color','r','FontUnits','normalized');        %scale letter to image

%8.rotate binary image by 45 degrees
I_bin_r=imrotate(I_bin,45,'crop');

%9.sum up a band of row (30%-60% of image size) to get row bndry for
row_band=zeros(1,size(I_bin_r,2));
for i=round(band_per(1,1)*size(I_bin_r,1)): round(band_per(1,2)*size(I_bin_r,1))
    row_band=row_band+double(I_bin_r(i,:));
end
ind_l_r=find(row_band>0,1);%v3
ind_r_r=find(row_band>0,1,'last');%v3

%10.sum up a band of col (30%-60% of image size) to get col bndry for
%  positive gradient length
col_band=zeros(size(I_bin_r,1),1);
for i=round(band_per(1,1)*size(I_bin_r,2)):...
        round(band_per(1,2)*size(I_bin_r,2))
    col_band=col_band+double(I_bin_r(:,i));
end
ind_t_r=find(col_band>0,1);%v3
ind_b_r=find(col_band>0,1,'last');%v3

%10.find pixel distance and convert to mm
distance_col_r=ind_r_r-ind_l_r;
distance_row_r=ind_b_r-ind_t_r;
distance_real_ng=distance_col_r*pxl_sz(2,1);
distance_real_ng=round(distance_real_ng*10)/10;
distance_real_pg=distance_row_r*pxl_sz(2,1);
distance_real_pg=round(distance_real_pg*10)/10;

ind_lr_r=(ind_b_r-ind_t_r)/2+ind_t_r;%v5
ind_tb_r=(ind_r_r-ind_l_r)/2+ind_l_r;%v5
subplot(1,2,2);
imshow(imrotate(I,45,'crop'),[]);
hold on
plot([ind_tb_r,ind_tb_r],[ind_t_r,ind_b_r],'Color','r','LineWidth',2);
plot([ind_l_r,ind_r_r],[ind_lr_r,ind_lr_r],'Color','r','LineWidth',2);
text(ind_lr_r-50,ind_t_r+20,...%HW:50 pxls to left
    [num2str(distance_real_ng) '\rightarrow'],...%use 'normalized' to
    'Color','r','FontUnits','normalized');       %scale letter to image
text(ind_lr_r+20,ind_tb_r-10,...%HW:20 pxls up
    ['\downarrow' num2str(distance_real_pg)],...%use 'normalized' to scale
    'Color','r','FontUnits','normalized');      %letter to image

end




