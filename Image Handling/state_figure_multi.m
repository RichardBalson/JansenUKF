function [fig_handle]= state_figure_multi(name,fig_structure,fig_settings,time,input,legendT,Error,Rows,Cols,yaxis,FigPlots)
% script created by Richard Balson 21/02/2013

EEG_Figure_multi;
fig_handle = figure('name',name,...
    'units','centimeters',...
    'position',[fig_settings.left_pos fig_settings.bottom_pos fig_width fig_height],...
    'papersize',[fig_width fig_height],...
    'filename',fig_dirandname,...
    'PaperPositionMode','auto');

for k = 1:size(input,1)/FigPlots
        subplot(Rows,Cols,k),plot(time,input((k-1)*FigPlots+1,:),color{1});
        hold on
if FigPlots>1
for j = 2:FigPlots
    plot(time,input((k-1)*FigPlots+j,:),color{j})
    hold on
end
if ~isempty(Error)
    plot(time,Error(k,:),ErrCol);
    hold on
    plot(time,Error(k+size(Error,1)/2,:),ErrCol);
end
end
ylabel(yaxis{k},'fontsize',fig_settings.label_fontsize)
set(gca,'fontsize',fig_settings.tick_fontsize)
box off
minc = min(min(input((k-1)*FigPlots+1:k*FigPlots,:)));maxc = max(max(input((k-1)*FigPlots+1:k*FigPlots,:)));
axis([0 max(time) (minc-abs(minc)*fig_settings.scale) (maxc+abs(maxc)*fig_settings.scale)]);
end
xlabel('Time (s)','fontsize',fig_settings.label_fontsize)
% title('Pyramidal Population Input','fontsize', label_fontsize)
k = legend(legendT,'Location',legLoc,'Orientation',legOri);
legend(k,'boxoff');
set(k,'fontsize',fig_settings.legend_fontsize);


