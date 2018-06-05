function [A,B,C,D,Distortion] = GlobalDistortion_S7(img,reso)

I = squeeze(img(:,:,6));
pxl_sz = reso(:);
per = .1;
diameter_in_mm = 190;
diamter_in_pixels = 190/pxl_sz(1);
radius_in_pixels = diamter_in_pixels/2;
conservative_radius_in_pixels = radius_in_pixels * 1.12;
[xp yp] = size(I);

I_max=double(max(max(I)));
[hist_cnt,hist_int]=hist(I(:),0:I_max);
hist_sample_start=round(I_max*per);%percentage of max as min to exclude air intenisty
%find max within specified histogram range
[int_cnt,int_pk]=max(hist_cnt(hist_sample_start:size(hist_cnt,2)));%v2
int_pk=int_pk+hist_sample_start-1;%v2
mu = int_pk;
% threshold at mu/2 of half water mean
I_bin = I > mu/2;
BW = imfill(I_bin,'holes');

[y, x] = ndgrid(1:size(BW, 1), 1:size(BW, 2));
centroid = mean([x(logical(BW)), y(logical(BW))]); % in x and y

% move from centroid in directions @30, 75, 120, 165 (and oppsoites) until
% reach .5 at BW_blend
figure; 
imshow(I,[]); hold on

phi = 30; phi_rad = phi*pi/180;
xi = centroid(1) + conservative_radius_in_pixels*cos(phi_rad);
yi = centroid(2) + conservative_radius_in_pixels * sin(phi_rad);
xf = centroid(1) + conservative_radius_in_pixels * cos(phi_rad + pi);
yf = centroid(2) + conservative_radius_in_pixels * sin(phi_rad + pi);
[cx,cy,c,~,~] = improfile(BW,[xi xf],[yi yf],1000);
% location of first and last point in mask
First = [cx(find(c>0,1)) cy(find(c>0,1))];
Last = [cx(find(c>0,1,'last')) cy(find(c>0,1,'last'))];
d_pixels = pdist([First;Last],'euclidean');
A = d_pixels*pxl_sz(1);
plot([First(1),Last(1)],[First(2),Last(2)],'Color','r','LineWidth',2);%v5

phi = 75; phi_rad = phi*pi/180;
xi = centroid(1) + conservative_radius_in_pixels*cos(phi_rad);
yi = centroid(2) + conservative_radius_in_pixels * sin(phi_rad);
xf = centroid(1) + conservative_radius_in_pixels * cos(phi_rad + pi);
yf = centroid(2) + conservative_radius_in_pixels * sin(phi_rad + pi);
[cx,cy,c,~,~] = improfile(BW,[xi xf],[yi yf],1000);
% location of first and last point in mask
First = [cx(find(c>0,1)) cy(find(c>0,1))];
Last = [cx(find(c>0,1,'last')) cy(find(c>0,1,'last'))];
d_pixels = pdist([First;Last],'euclidean');
B = d_pixels*pxl_sz(1);
plot([First(1),Last(1)],[First(2),Last(2)],'Color','r','LineWidth',2);%v5

phi = 120; phi_rad = phi*pi/180;
xi = centroid(1) + conservative_radius_in_pixels*cos(phi_rad);
yi = centroid(2) + conservative_radius_in_pixels * sin(phi_rad);
xf = centroid(1) + conservative_radius_in_pixels * cos(phi_rad + pi);
yf = centroid(2) + conservative_radius_in_pixels * sin(phi_rad + pi);
[cx,cy,c,~,~] = improfile(BW,[xi xf],[yi yf],1000);
% location of first and last point in mask
First = [cx(find(c>0,1)) cy(find(c>0,1))];
Last = [cx(find(c>0,1,'last')) cy(find(c>0,1,'last'))];
d_pixels = pdist([First;Last],'euclidean');
C = d_pixels*pxl_sz(1);
plot([First(1),Last(1)],[First(2),Last(2)],'Color','r','LineWidth',2);%v5

phi = 165; phi_rad = phi*pi/180;
xi = centroid(1) + conservative_radius_in_pixels*cos(phi_rad);
yi = centroid(2) + conservative_radius_in_pixels * sin(phi_rad);
xf = centroid(1) + conservative_radius_in_pixels * cos(phi_rad + pi);
yf = centroid(2) + conservative_radius_in_pixels * sin(phi_rad + pi);
[cx,cy,c,~,~] = improfile(BW,[xi xf],[yi yf],1000);
% location of first and last point in mask
First = [cx(find(c>0,1)) cy(find(c>0,1))];
Last = [cx(find(c>0,1,'last')) cy(find(c>0,1,'last'))];
d_pixels = pdist([First;Last],'euclidean');
D = d_pixels*pxl_sz(1);
plot([First(1),Last(1)],[First(2),Last(2)],'Color','r','LineWidth',2);%v5

distances = [A B C D];
Distortion = 100*(max(distances)-min(distances))/mean(distances);






