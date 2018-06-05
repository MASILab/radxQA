function [] = print_current_figure(rez,name)
% prints current figure @ given resolution and with given name

% rez=400; %resolution (dpi) of final graphic

f=gcf; %f is the handle of the figure you want to export

figpos=getpixelposition(f); 

resolution=get(0,'ScreenPixelsPerInch'); 

set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 

print(f,name,'-dpng',['-r',num2str(rez)],'-opengl') %save file 

end
