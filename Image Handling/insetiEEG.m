function [h_main, h_inset]=insetiEEG(main_handle, inset_handle,inset_size)

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.

if nargin==2
    inset_size=0.35;
end

inset_size=inset_size*.7;
% figure
figure('color','white','units','centimeters','position',[2 2 8.6 3],'papersize',[8.6 3],'filename','iEEG.pdf')
new_fig=gcf;
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);
set(h_main,'Position',get(main_fig,'Position'),'box','off')
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax=get(main_fig,'Position');
set(h_inset,'Position', [.7*ax(1)+ax(3)-inset_size .5*ax(2)+ax(4)-inset_size+.1 inset_size-.1 inset_size])
FS_labels=7;
leg=legend(h_main,'$v_p(t)$')
set(leg,'interpreter','latex')
set(leg,'fontsize',FS_labels,'box','off','Orientation','horizontal', 'units','centimeters', 'position', [1 2.6 6.6 0.2],'LineWidth', .1) 

