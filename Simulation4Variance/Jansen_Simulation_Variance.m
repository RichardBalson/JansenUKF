function Jansen_Simulation_Variance

% script created by Richard Balson 21/02/2013

% description
% ~~~~~~~~~~~
% Simulate states and output of the extended neural mass model for
% estimation purposes. Input mean assumed to be 90.

% last edit
% ~~~~~~~~~
clear
close all
clc

% Variables required
% ~~~~~~~~~~~~~~~
Simulation_number =1000;
addpath(genpath('../../Jansen')); % Specify files required for estimation
% stochastic -  Determines level of stochasticity in input
% number_of_sigma_input - Determines the standard deviation of stochastic input
% Parameter_index - Used to decide what parameter should be simulated
% Input_mean_variation - uSed to specifiy whether the mean of the input
% should be different from its described nominal value
SimulationSettings.name = 'VarianceSim'; % Specify name of file to save Wendling model output data to, or to load data from when simulation is not performed
SimulationSettings.fs = 2048;
SimulationSettings.simulation_time =40; %Time for simulation in seconds 
    SimulationSettings.slope_time =2.5; % Specifies the time over which the model gain should be altered
    SimulationSettings.number_of_sigma_input = 4; % Used to determine standard deviation of input if  1: 68.27% of realisations within physiolgical range, 2: 95.45, 3: 99.73 4: 99.994
    SimulationSettings.stochastic = 1; % Used to specifiy the stochastic adjustment on the input 1 is no adjustment. <1 downscalling, >1 upscaling
    SimulationSettings.Parameter_index = 6; % Choose parameters to be simulated: 1 = Seizure Parameter from Wendling 2002;
    %  2 = Seizure Parameter from Wendling 2005;...
    %  3 = Altered excitability;
    %  4 = Parameters at midpoint of their range;
    %  5 = random parameters;
    %  6= random parameters and random number of
    %  variations in simulation
    % 7 User defined
    if SimulationSettings.Parameter_index % Specify synaptic gains for simulation
        SimulationSettings.AV= [4 5 6 7 8];
        SimulationSettings.BV= [20 15 10 5 0];
    end
    SimulationSettings.Input_mean_variation = 0; % If 0 mean stays constant for simulation,
    %if 1 input mean is drawn from a uniform distribution limited by the physiological limits of the input
    % if 2 input mean is drawn from a Gaussian
    % distribution with a mean as per Wendling
    % 2002 and stnadard deviation that satisfies
    % number_of_sigma_input
    
% next edit
% ~~~~~~~~~

% Beginning of script
% ~~~~~~~~~~~~~~~~~~~~~

% Model parameters(Variable parameters are to alter parameter values during
% simulation) Matrix specifies parameter values for a specified time period
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Units mV
% ~~~~~~~~~~~~

sampling_frequency = SimulationSettings.fs; % Sampling frequency for solutions and measurements
Max_parameters = SimulationSettings.simulation_time/(SimulationSettings.slope_time+0.5); % Specifiy maximum number of parameters for random simulation

for count=1:Simulation_number
[AV BV] = Parameter_choice(SimulationSettings,Max_parameters);

% Solver parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simulation_changes = length(AV); % Specify the number of times the parameters change in simulation

% Units Hz
% ~~~~~~~~


%%
% Stationery Parameter defintions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~
JansenModelStationaryParameters


% e = xxxx*randn(NSamples,1);
% Determine stochastic input to the model
% ~~~~~~~~~~~~~~~~~~~~~~~
gauss = randn(1,length(normalise)+1)*std_deviation*SimulationSettings.stochastic;

%         gauss = randn(1)*std_deviation*stochastic; % Determine random fluctuations in input signal with the specified standard deviation

normalised_gaussian_input(1,:) = meanf + gauss(1,:); % Determine value of input for all time steps

%         while (any(normalised_gaussian_input>max_frequency) || any(normalised_gaussian_input <min_frequency))
%
%            gauss(1,k) = randn(1)*std_deviation*stochastic;
%
%             normalised_gaussian_input(k) = meanf + gauss(1,k);
%
%         end

normalised_gaussian_inputSDE(1,:) = meanf + gauss(1,:)*sqrt(dt)*sampling_frequency; % Determine value of input used for calculation purposes for all time steps

%%

% Equations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Solver
% ~~~~~~~~~~~~~~~~~~~~~~~
% Iteration to change parameters in function
for k = normalise
    %
    % 8th Order
    % Derivation of equation shown in pdf attached
    
    % Normalised gaussian input
    
    zs1 =  z(k,2)-z(k,3);
    zs2 = C(1)*z(k,1);
    zs3 = C(3)*z(k,1);
    [z(k+1,1) z(k+1,4)] = PSPkernel([z(k,1); z(k,4)],dt,MVI(k,1),tcon(1),sigmoid(zs1));
    [z(k+1,2) z(k+1,5)] = PSPkernel([z(k,2); z(k,5)],dt,MVI(k,1),tcon(1),normalised_gaussian_inputSDE(k) +C(2)*sigmoid(zs2));
    [z(k+1,3) z(k+1,6)] = PSPkernel([z(k,3); z(k,6)],dt,MVI(k,2),tcon(2),C(4)*sigmoid(zs3));
    
end% End of for loop

% Determine output of model
% ~~~~~~~~~~~~~~~~~~~~~~~~
output6 = z(:,2)-z(:,3);      % 8th order output, Normalised gaussian input

meanX(:,count) = mean(z);
covX(:,:,count) = cov(z);
stdX(:,count) = std(z);
end

save([SimulationSettings.name,'.mat'],'meanX', 'stdX','covX');
