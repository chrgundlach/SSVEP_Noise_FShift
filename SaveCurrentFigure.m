function [  ] = SaveCurrentFigure( path, name )
%SaveCurrentFigure Saves current figure with standard values
%   saves figure with 600 dpi
%   as .fig, .eps, .jpeg


%% check inputs
if ~isdir(path)
    mkdir(path)
    fprintf(1,'\n...created directory: %s',path)
end

%% set image properties
set(gcf,'PaperPositionMode','auto','Color',[1 1 1])

% set font
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')

try set(gca,'fontname','Arial'); end
try set(gcf,'fontname','Arial'); end
h_child = get(gcf,'children');
for i_ch=1:length(h_child)
    h_text = findobj(get(h_child(i_ch),'children'),'-not','FontName','Arial');
    try set(h_child(i_ch),'FontName','Arial'); end
    try set(h_text,'FontName','Arial'); end
end

%% save image
fprintf(1,'\n...saving files in path: %s',path)
% as jpeg
fprintf(1,'\n...%s.png',name)
print(fullfile(path,name),'-dpng','-r600')
% as fig
fprintf(1,'\n...%s.fig',name)
saveas(gcf,fullfile(path,name),'fig')
% as eps
fprintf(1,'\n...%s.eps',name)
% print(fullfile(path,name),'-depsc','-r600')
saveas(gcf,fullfile(path,name),'epsc')
% exportgraphics(gcf,fullfile(path,name),'BackgroundColor','none','ContentType','vector')
fprintf(1,'\n...done!\n')

end

