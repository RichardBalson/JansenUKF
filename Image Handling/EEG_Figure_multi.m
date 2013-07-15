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

if strcmp(fig_structure,'State')
fig_width = 13*Cols;                % cms
fig_height = 4*Rows;                 % cms
fig_dirandname = [fig_settings.dirname name];
legOri = 'horizontal';
color = {'k' 'r' 'b' 'g'};
ErrCol = 'g';
legLoc = [0.3 0.45 0.5 1];
linewidth =1;
elseif strcmp(fig_structure,'Multi')
fig_width = 13*Cols;                % cms
fig_height = 4*Rows;                 % cms
fig_dirandname = [fig_settings.dirname name];
legLoc = [0.3 0.45 0.5 1];
legOri = 'horizontal';
color = {'k' 'r' 'b' 'g'};
ErrCol = 'g';
linewidth =2.5;
end