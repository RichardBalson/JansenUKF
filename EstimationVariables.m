% script created by Richard Balson 21/02/2013

% description
% ~~~~~~~~~~~
% this script assignes all static variables, all staes and parameters are
% initialised by realising a gaussin distribution with an assigned number
% of standard deviations.

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~

% Beginning of script

% ~~~~~~~~~~~~~~~~~~~~~

% Image names and handling
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~`
if ((min(MVI(:,1))==max(MVI(:,1))) && (min(MVI(:,2))==max(MVI(:,2))))
    Parameter_varying=0; % Initialse parameter m, If m = 1 parameter uncertainty increases to account for the fact that the parameters vary in time. m is always set to zero, and is adjusted in the code accordingly
    simulation_initial_name = [Estimation_Type,'f',int2str(SimulationSettings.fs),...
        'A,B',int2str(MVI(1,1)),',',int2str(MVI(1,2)),'S',num2str(SimulationSettings.stochastic)]; % Initaite name for simulation, used for saving purposes
else
    Parameter_varying=1;
    simulation_initial_name = [Estimation_Type,'Vf',int2str(SimulationSettings.fs),...
          'A,B',int2str(MVI(1,1)),',',int2str(MVI(1,2)),'S',num2str(SimulationSettings.stochastic)];  % Initaite name for simulation, used for saving purposes
end

sampling_frequency =SimulationSettings.fs;
% Assign standard deviation for all variables
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% load InitialP1E CheckPxx CheckMean CheckInit
% 
% mean_state = (mean(CheckInit,2))'; % Deterine mean based on numerous previous simulation means
% 
% State_std_deviation = (sum(CheckPxx(1:8,1:8,:),3)/(size(CheckPxx,3)-1)*ones(8,1))';
% 
% Base_state_uncertainty = (State_std_deviation.^2)*sqrt(dt); % Determine base uncertainty for states

load VarianceSim meanX stdX


mean_state = (mean(meanX,2))'; % Deterine mean based on numerous previous simulation means

State_std_deviation = (sum(stdX,2)/(size(stdX,2)-1))';

Base_state_uncertainty = (State_std_deviation.^2)*sqrt(dt); % Determine base uncertainty for states


% Determine standard deviation for all slow states
% ~~~~~~~~~~~~~~~~~~~~~~~~~~

State_sigma = State_std_deviation/number_of_sigma;

% Model Input
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Input_mean = (frequency_limits(2)+frequency_limits(1))/2; % Mean input frequency from pyramidal neruons external to model

% Simulated Signal Data mV
% ~~~~~~~~~~~~~~~~~~
EstStart_Sample = floor(EstStart*sampling_frequency);
check = output6(EstStart_Sample:end); % Assign a new variable as the simulated output, data point form time EstStart are used as the beginning of the observations
Y = check + randn(length(check),1)*NoiseIn; % Initalise noise on the simulated signal
Number_of_observations = length(Y); % Parameter that defines he number of observations from the simulation

% Number of sigma points
% ~~~~~~~~~~~~~~~~~~~

if kappa > 0
    Sigma_points = 2*Dx+1;  % Specify number of sigma points to generate, if kappa is greater than zero than the mean is propagated
    % through the system and there is therefore one more sigma point generated
else
    Sigma_points = 2*Dx; % When kappa is set to zero sigma points are generated, note there are twice the number of sigma points as there are states
    % the reason for this is that for each mean state
    % value, a sigma point one standard deviation from
    % its mean are propagated through the system. This
    % inclue the sigma points: mean - standard deviation
    % and mean + standard deviation
end

% Estimation states uncertainty
% ~~~~~~~~~~~~~~~~~~~~~~~~~

% Uncertianty for each fast state
% mV

State_uncertainty = Base_state_uncertainty./State_uncertainty_adjustment; % Specify base
% state uncertainty for all states
% State_uncertainty(1,[2 6]) = (ones(1,2)*stochastic*Variable_state_uncertainty+State_uncertainty(1,[2 6])); % Alter uncertainty of parameter affected directly by stochastic input

% Uncertainty for slow states
% ~~~~~~~~~~~~~~~~~~~~~~~

if Dp >0
    Parameter_uncertainty(1,1) = Exc_parameter_uncertainty + Variable_parameter_uncertainty*Parameter_varying; % Variance to allow for model error and stochastic input effect
    if Dp >1
        Parameter_uncertainty(1,2) = (Base_parameter_uncertainty + Variable_parameter_uncertainty*Parameter_varying);
    end
else
    Parameter_uncertainty = [];
end

% Uncertainty for model input mean
% ~~~~~~~~~~~~~~~~~

if Dk==1
    Input_uncertainty = Variable_input_uncertainty*(SimulationSettings.Input_mean_variation>0)+Base_input_uncertainty; % Define the uncertainty of the input to the model, 
    Input_uncertainty = Input_uncertainty.^2.*uncertainty_adjustment;
else
    Input_uncertainty =[];
end

% Matrix manipulation
% ~~~~~~~~~~~~~~~~~~~~~~~~

State_sigma = repmat(State_sigma,Ds,1);

State_uncertaintyF = repmat(State_uncertainty,Ds,1);

Parameter_uncertaintyF = repmat(Parameter_uncertainty,Dp,1);% Define a parameter uncertainty matrix that can easily be used for future purposes.


%%%%%%%%%%%%%%%%%%%%%%%%%
% Intialise all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(Dx,Number_of_observations); % Intialise state estimate matrix
Pxx = zeros(Dx,Dx,Number_of_observations);% Intialise state covariance estimate matrix
Xout = zeros(Dx,Sigma_points,Number_of_observations); % Intialise output states from the point sin the unscented transform matrix
Yout = zeros(Dy,Sigma_points,Number_of_observations); % Intialsie the model output for each set of sigm points matrix
ExpX = zeros(Dx,Number_of_observations);% Intialise athe expected states from the unscented transform matrix
ExpY = zeros(Dy,Number_of_observations);% Intialise the expected states from the unscented transform matrix
Sigma = zeros(Dx,Sigma_points,Number_of_observations);% Intialise the sigma points matrix

State_covariance_matrix = eye(Ds).*State_sigma.^2; % State covariance matrix

State_uncertainty_matrix = eye(Ds).*State_uncertaintyF.^1.75.*uncertainty_adjustment;

Parameter_uncertainty_matrix = eye(Dp).*Parameter_uncertaintyF.^2.*uncertainty_adjustment;

R = Observation_uncertainty^2;

Q = blkdiag(State_uncertainty_matrix,Input_uncertainty,Parameter_uncertainty_matrix);
