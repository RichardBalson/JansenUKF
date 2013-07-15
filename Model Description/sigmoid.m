% script created by Richard Balson 05/10/2011

% description
% ~~~~~~~~~~~
% This function describes a sigmoid function which convertes membrane
% potentials into neuronal spike trains. The input to the model is the the
% specified membrane potential in mV. The output is the number of spikes
% per second.

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function s = sigmoid(v)

% Parameter specifiction
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Units seconds^(-1)
% ~~~~~~~~~~~

e =2.5;      %sigmoid parameters, maximal firing rate of neuronal population

% Units mV
% ~~~~~~~~~~~~

k =6;        %v0  PSP at which 50% of maximum firing rate is achieved

% Units mV^(-1)
% ~~~~~~~~~~~

r =.56; % Indicates the steepness of the sigmoid funtion
    
% Equations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

s = 2*e./(1+exp(r*(k-v))); % Equation to convert PSP into a firing rate


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% End of function description

