function varargout=IMcircle(sizes,center,radius)
%_______________________________________________________________________
%function IMcircle(sizes,center,radius);
%_______________________________________________________________________
%
% This function creates a 2D binary image of a circle or group of circles.
% The mandatory inputs are the size of the 2D image output, the coordinates
% of the center of the circle (i.e. [x,y]), and the radius of the circle in
% units of pixels.  Each of these pieces of information should be organized
% along the columns of the inputs.  Multiple circles can be drawn by adding
% new entries along the rows of these inputs.  The output is a 2D matrix of
% size 'sizes' with values of zero everywhere outside the circle/circles,
% and ones inside the circle/circles.
%
% Example: imagesc(IMcircle([500,500],[100,130],50));
%          axis equal;
% Explanation: This will create a 500x500 pixel image with zeros everywhere
%           except a circle centered at [100,130] whose radius is 50 pixels
%           wide.
%
% Written by: Allen T. Newton, PhD   20120419
%_______________________________________________________________________
%

%---------------------
%check inputs
%---------------------

if nargin~=3
    error('Required inputs ''sizes'',''center'', and ''radius'' must be all supplied');
end
if size(sizes,2)~=2
    error('''sizes'' must be a two element column-wise vector');
end
if size(center,2)~=2
    error('''center'' must be a two element column-wise vector');
end
if ~isequal(size(center,1),size(radius,1))
    error('inputs ''sizes'', and ''radius'' must be the same length along the first dimension');
end


obj=zeros(sizes);
%loop through each specified circle
for ii=1:size(center,1);
    tmp=zeros(size(obj));
    tmp(center(ii,1),center(ii,2))=1; 
    circ_mask = double(getnhood(strel('ball',radius(ii),radius(ii),0)));
    tmp=conv2(tmp,circ_mask,'same');
    obj=obj+tmp;
end
clear tmp;




obj=logical(obj);
%obj=obj(1:sizes(1),1:sizes(2))~=0;
varargout{1}=obj;
% imagesc(obj);axis image;colormap gray;hold all;plot(x,y1,'r',x,y2,'g');hold off;





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%   PREVIOUS VERSION    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% 
% function varargout=IMcircle(sizes,center,radius)
% %_______________________________________________________________________
% %function IMcircle(sizes,center,radius);
% %_______________________________________________________________________
% %
% % This function creates a 2D binary image of a circle or group of circles.
% % The mandatory inputs are the size of the 2D image output, the coordinates
% % of the center of the circle (i.e. [x,y]), and the radius of the circle in
% % units of pixels.  Each of these pieces of information should be organized
% % along the columns of the inputs.  Multiple circles can be drawn by adding
% % new entries along the rows of these inputs.  The output is a 2D matrix of
% % size 'sizes' with values of zero everywhere outside the circle/circles,
% % and ones inside the circle/circles.
% %
% % Example: imagesc(IMcircle([500,500],[100,130],50));
% %          axis equal;
% % Explanation: This will create a 500x500 pixel image with zeros everywhere
% %           except a circle centered at [100,130] whose radius is 50 pixels
% %           wide.
% %
% % Written by: Allen T. Newton, PhD   20120419
% %_______________________________________________________________________
% %
% 
% %---------------------
% %check inputs
% %---------------------
% 
% if nargin~=3
%     error('Required inputs ''sizes'',''center'', and ''radius'' must be all supplied');
% end
% if size(sizes,2)~=2
%     error('''sizes'' must be a two element column-wise vector');
% end
% if size(center,2)~=2
%     error('''center'' must be a two element column-wise vector');
% end
% if ~isequal(size(center,1),size(radius,1))
%     error('inputs ''sizes'', and ''radius'' must be the same length along the first dimension');
% end
% 
% 
% obj=zeros(sizes);
% %loop through each specified circle
% for ii=1:size(center,1);
% 
%     x=[1:sizes(2)];
%     
%     y1=(floor(real(center(ii,2)-sqrt((radius(ii).^2)-((x-center(ii,1)).^2)))))+1;
%     y2=(ceil(real(center(ii,2)+sqrt((radius(ii).^2)-((x-center(ii,1)).^2)))))-1;
% 
%     %keyboard;
%     for jj=1:size(obj,2);
%         obj(y1(jj):y2(jj),jj)=1;
%         obj(y2(jj):y1(jj),jj)=0;
% %         obj(center(ii,2),x(find(y2~=y1,1,'first')-1))=1;
% %         obj(center(ii,2),x(find(y2~=y1,1,'last')+1))=1;
%     end    
% end
% obj=logical(obj);
% %obj=obj(1:sizes(1),1:sizes(2))~=0;
% varargout{1}=obj;
% % imagesc(obj);axis image;colormap gray;hold all;plot(x,y1,'r',x,y2,'g');hold off;