% script created by Richard Balson 29/03/2013

% description
% ~~~~~~~~~~~
% This script assigns all intial values for all model states

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~

% Beginning of script

% ~~~~~~~~~~~~~~~~~~~~~


% Initialise model slow states
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (Parameter_initialisation ==0)
    if Random_number_generator(2) % random number generator for intialisation of parameters
    rng(0);
    else
    rng(cputime);
    end
    mu_A = (Max_A+Min_A)/2;
    mu_B = (Max_B+Min_B)/2;
    mu_Input = (frequency_limits(2)+frequency_limits(1))/2;
    
    % Determine one standard deviation for all parameters
    
    sigma_A = (Max_A-Min_A)/(2*number_of_sigma);
    sigma_B = (Max_B-Min_B)/(2*number_of_sigma);
    sigma_Input = (frequency_limits(2)-frequency_limits(1))/(2*number_of_sigma);
        conditionP=1;
    while conditionP
    Initial_excitability = mu_A + randn(1)*sigma_A; % Intial excitability, estimate of parameter A
    Initial_inhibition = mu_B + randn(1)*sigma_B; % Intial slow inhibition, estimate of parameter B
    Initial_Input = mu_Input + randn(1)*sigma_Input;
       conditionP = any([Initial_excitability Initial_inhibition Initial_Input]<0);
    end
else
    Initial_excitability = (1-PercError/100)*MVI(1,1); % Intial excitability, estimate of parameter A
    Initial_inhibition = (1-PercError/100)*MVI(1,2); % Intial slow inhibition, estimate of parameter B
    Initial_Input = (1-PercError/100)*Input_mean;
end

% Intialise model fast states

Gauss_state = mean_state + State_sigma(1,:).*randn(1,Ds);

% Determine standard deviatio for all slow states
% ~~~~~~~~~~~~~~~~~~~~~~~~~

if Dp >0
    Parameter_std_deviation(1,1) = max(Max_A-Initial_excitability,Initial_excitability-Min_A)/number_of_sigma*std_adjustment_Exc; % Variance to allow for tracking of parameters Parameters appear to work  5e-3
    if Dp >1
        Parameter_std_deviation(1,2) = max(Max_B-Initial_inhibition,Initial_inhibition-Min_B)/number_of_sigma*std_adjustment_SInh;
    end
else
    Parameter_std_deviation = [];
end

% Initialise model input mean as a slow state
% ~~~~~~~~~~~~~~~~~~~

% Uncertainty for model input mean
% ~~~~~~~~~~~~~~~~~

if Dk==1
    X(Ds+1) = Initial_Input; % Initialise guess for input mean
    Input_std_deviation = max((X(Ds+1) - frequency_limits(1)),frequency_limits(2)-X(Ds+1))/std_adjustment_input; % Define the variance of the input to the model Base 1e-2
    Input_variance = Input_std_deviation.^2;
else
    Input_variance =[];
end


% Matrix Manipulation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~

Parameter_std_deviationF = repmat(Parameter_std_deviation,Dp,1);% Define a parameter variance matrix that can easily be used for future purposes.

% Initialise state and covariance matrix
% ~~~~~~~~~~~~~~~~~~~

X(1:Ds,1) = Gauss_state; % Intialise all parameters for the state estimate

if Dp >0                                    % Check whether parameters need to be estimated, if so alter storage folder
    Estimation_type = 'Results\ParameterEstimation'; % Set folder into which images will be sved if they are printed
    X(Ds+Dk+1,1) = Initial_excitability;     % Set the intial guess for the A to be the initial value specified.
    if Dp >1                                % Check whether more than one parameter needs to be estimated
        X(Ds+Dk+2,1) = Initial_inhibition;    % Set the intial guess for the B to be the initial value specified.
    end
else
    Estimation_type = 'Results\StateEstimation'; %Set folder into which images will be sved if they are printed
end

Parameter_covariance_matrix = eye(Dp).*Parameter_std_deviationF.^2; % State Parameter matrix

Pxx(:,:,1) = blkdiag(State_covariance_matrix,Input_variance,Parameter_covariance_matrix); %Covariance of model states and parameter

