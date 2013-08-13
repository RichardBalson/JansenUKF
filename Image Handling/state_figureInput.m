function [fig_handle]= state_figureInput(name,fig_structure,fig_settings,time,input,legendT,Error)
% script created by Richard Balson 21/02/2013

EEG_FigureF;
fig_handle = figure('name',name,...
    'units','centimeters',...
    'position',[fig_settings.left_pos fig_settings.bottom_pos fig_width fig_height],...
    'papersize',[fig_width fig_height],...
    'filename',fig_dirandname,...
    'PaperPositionMode','auto');


h=plot(time,input(1,:),color{1});
hold on
for k = 2:size(input,1)
    plot(time,input(k,:),color{k})
    hold on
end
if ~isempty(Error)
    m=plot(time,Error(1,:),ErrCol);
    hold on
    n=plot(time,Error(2,:),ErrCol);
    
end

set(gca,'fontsize',fig_settings.tick_fontsize)
box off
xlabel('Time (s)','fontsize',fig_settings.label_fontsize)
ylabel('Frequency (Hz)','fontsize',fig_settings.label_fontsize)
% title('Pyramidal Population Input','fontsize', label_fontsize)
k = legend(legendT,'Location',legLoc,'Orientation',legOri);
legend(k,'boxoff');
uistack(h,'top');
if ~isempty(Error)
    uistack(n,'bottom');
    uistack(m,'bottom');
end
set(k,'fontsize',fig_settings.legend_fontsize);
minc = min(min(input));maxc = max(max(input));
axis([min(time) max(time) (minc-abs(minc)*fig_settings.scale) (maxc+abs(maxc)*fig_settings.scale)]);

