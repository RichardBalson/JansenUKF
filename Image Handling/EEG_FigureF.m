% script created by Richard Balson 22/02/2013

% description
% ~~~~~~~~~~~

% This script specifies the setting for an EEG image

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~

% Beginning of script
% ~~~~~~~~~~~~~~~~~~~~~

if strcmp(fig_structure,'Obs')

fig_width = 14.5;                % cms
fig_height = 5;                 % cms
fig_dirandname = [fig_settings.dirname name];
legLoc= 'NorthOutside';
legOri = 'horizontal';
color = {'b' 'r' 'k' 'g'};

elseif strcmp(fig_structure,'State')
fig_width = 14.5;                % cms
fig_height = 5;                 % cms
fig_dirandname = [fig_settings.dirname name];
legLoc= 'NorthOutside';
legOri = 'horizontal';
color = {'k' 'r' 'b' 'g'};
ErrCol = 'g';
end