% script created by Richard Balson 06/09/2012

% description
% ~~~~~~~~~~~
% This function describes the Kalman filter

% last edit
% ~~~~~~~~~

% Calculation of covariacne matrix altered from Pxxf = Pxx-K*Pxy' to
% K*Pyy*K'

% next edit
% ~~~~~~~~~

% Beginning of function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [X Pxxf] = Kalman(ExpX, ExpY, Ym, Pxx, Pxy, Pyy)

K = Pxy/Pyy;
X = ExpX + K*(Ym-ExpY);
Pxxf = Pxx - K*Pxy'; % K*Pyy*K' or K*Pxy'
 
end

% Parameter specifiction
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Units
% ~~~~~~~~

% Equations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% End of function description
