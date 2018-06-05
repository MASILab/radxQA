function [centre_coord,contour_xy,I_bin]=MakeEllipseROI(I,phan_centre,phan_bndry,FOV_bndry,wl_pxl,pxl_sz,choice)

%1.check if there is enough space to create ROI between phantom and FOV
empty_space=abs(phan_bndry-double(FOV_bndry));
ROI_phan_distance=5;%HW:default ROI to phantom distance
if empty_space>wl_pxl(1,1)+ROI_phan_distance+1%v2
    disp('You have enough space to create ROI.');
    w_pxl=wl_pxl(1,1);
    l_pxl=wl_pxl(1,2);
else
    disp('Not enough empty space, re-calculate ROI width to fit in.');
    w_pxl=(empty_space-ROI_phan_distance-1);%1 pxl away from FOV edge%v2
    l_pxl=4*1000/(pi*w_pxl)/pxl_sz(2,1);%make sure 10cm2 area%v2
    fprintf('New ROI width is %2d and ROI length is %2d\n',w_pxl,l_pxl);
end
%2.define ellipse ROI centre coord
switch choice
    case 2
        x_c=phan_centre(1,2);
        y_c=phan_bndry-ROI_phan_distance-w_pxl/2+1;%1 pxl away from FOV edge
    case 3
        x_c=phan_centre(1,2);
        y_c=phan_bndry+ROI_phan_distance+w_pxl/2+1;%1 pxl away from FOV edge
    case 4
        x_c=phan_bndry-ROI_phan_distance-w_pxl/2+1;%1 pxl away from FOV edge
        y_c=phan_centre(1,1);
    case 5
        x_c=phan_bndry+ROI_phan_distance+w_pxl/2+1;%1 pxl away from FOV edge
        y_c=phan_centre(1,1);
end
centre_coord=[x_c,y_c];
%3.create mask image
theta=0:0.01:2*pi;
if choice==2 || choice==3
    x=l_pxl/2*cos(theta)+x_c;%v2
    y=w_pxl/2*sin(theta)+y_c;%v2
elseif choice==4 || choice==5
    x=w_pxl/2*cos(theta)+x_c;%v2
    y=l_pxl/2*sin(theta)+y_c;%v2
end
I_bin=roipoly(I,x,y);
%4.output x&y contour
contour_xy=[x;y];
end