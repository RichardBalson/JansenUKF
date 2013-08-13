function [fig_handle]= state_figure_multi_sim(name,fig_structure,fig_settings,time,input,legendT,Error,ErrorP,Rows,Cols,yaxis,FigPlots)
% script created by Richard Balson 21/02/2013
% Check that input is a column vector...................

EEG_Figure_multi;
fig_handle = figure('name',name,...
    'units','centimeters',...
    'position',[fig_settings.left_pos fig_settings.bottom_pos fig_width fig_height],...
    'papersize',[fig_width fig_height],...
    'filename',fig_dirandname,...
    'PaperPositionMode','auto');

for k = 1:size(input,3)
    subplot(Rows,Cols,k),h=plot(time,input(1,:,k),color{1},'LineWidth',linewidth);
    hold on
    if size(input,1)>1
        plot(time,input(2:end,:,k))
        hold on
        if ~isempty(Error)
            m(k)=plot(time,Error(1:end,:,k),ErrCol);
            hold on
            n(k)=plot(time,ErrorP(1:end,:,k),ErrCol);
            
        end
    end
    ylabel(yaxis{k},'fontsize',fig_settings.label_fontsize)
    set(gca,'fontsize',fig_settings.tick_fontsize)
    box off
    minc = min(min(input(:,:,k)));maxc = max(max(input(:,:,k)));
    axis([min(time) max(time) (minc-abs(minc)*fig_settings.scale) (maxc+abs(maxc)*fig_settings.scale)]);
    
end
xlabel('Time (s)','fontsize',fig_settings.label_fontsize)
% title('Pyramidal Population Input','fontsize', label_fontsize)
k = legend(legendT,'Location',legLoc,'Orientation',legOri);
if ~isempty(Error)
    for a =1:size(input,3)
        uistack(n(a),'bottom');
        uistack(m(a),'bottom');
    end
end
uistack(h,'top');
legend(k,'boxoff');
set(k,'fontsize',fig_settings.legend_fontsize);


