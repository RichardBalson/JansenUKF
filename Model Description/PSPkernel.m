% script created by Richard Balson 05/10/2011

% description
% ~~~~~~~~~~~
% This function describes the equation used for the PSP kernel, this
% function converts firing rates into PSP

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~

% ~~~~~~~~~~~~~~~~~~~~~~~


function [x y]= PSPkernel(yin,dt,Gain,tconstant,input) % Function to convert firing rate to PSD for each equation, yin is the current solutions input,
% dt the period of sampling, Gain the gain parameter used for the PSP,
% tconstant the time constant used by the PSP, input is the input to the
% system. x and y indicate the output of the model
                                                       % 
%%

x = yin(1,:) + yin(2,:)*dt; % dx/dt = y
y = yin(2,:) + (Gain*tconstant.*input - 2*tconstant*yin(2,:) -tconstant^2*yin(1,:))*dt; % dy/dt = A*a*input-2*a*y-a^2*x; where tconstant and gain have been replaced by a and A respectively

end



