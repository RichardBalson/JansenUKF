% script created by Richard Balson 09/01/2013

% description
% ~~~~~~~~~~~
% this script determines numerous error results from UKF estimation of
% states and parameters

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~

% Beginning of script
% ~~~~~~~~~~~~~~~~~~~~~

% Determine mean square error for states
% ~~~~~~~~~~~~~~~~~~

limit = length(Y);

errMulti = zeros(Dp+Dk,Simulation_number); % Initlaise mean square error
NerrMulti = zeros(Dp+Dk,Simulation_number); % Initialise normalised error for initial sum
PercErrEstimateMulti =zeros(Dp+Dk,size(X_Multi,1),Simulation_number);

    Sample_inverse = 1/limit; % Used to normalise mean squared error
    
    IndexErrD = (EstStart)*sampling_frequency+1:Decimate:length(output6);

for q = 1:Simulation_number
    
    % Determine mean square error for parameter estimates
    
    if (Dp)>0
        errMulti(1:Dp,q) = diag(Sample_inverse*(MVI(IndexErrD,1:Dp) - X_Multi(:,Dk+1:Dk+Dp,q))'*(MVI(IndexErrD,1:Dp) - X_Multi(:,Dk+1:Dk+Dp,q)));
    end
    
    if Dk >0
        errMulti(Dp+1,q) = Sample_inverse*(ones(length(IndexErrD),1)*Input_mean-X_Multi(:,1,q))'*(ones(length(IndexErrD),1)*Input_mean - X_Multi(:,1,q));
    end
    
    % Determine normalised error
    
    if (Dp)>0
        New = MVI(IndexErrD,1:Dp);
        New(New==0) = eps;
        NerrMulti(1:Dp,q) = Sample_inverse*diag((MVI(IndexErrD,1:Dp) - X_Multi(:,Dk+1:Dk+Dp,q))'./New(:,1:Dp)'*(MVI(IndexErrD,1:Dp) - X_Multi(:,Dk+1:Dk+Dp,q)));
        PercErrEstimateMulti(1:Dp,:,q) = (abs(X_Multi(:,Dk+1:Dk+Dp,q)-MVI(IndexErrD,1:Dp))./MVI(IndexErrD,1:Dp)*100)';
    end
    
    if Dk >0
        NerrMulti(Dp+1,q) = Sample_inverse*diag((Input_mean-X_Multi(Ds+1,:))/Input_mean*(Input_mean-X_Multi(Ds+1,:))');
        PercErrEstimateMulti(Dp+1,:,q) = (abs(X_Multi(:,1,q)-Input_mean)./Input_mean*100)';
    end
end

Error_name = strcat(Estimation_type,'\Error_Multi_',simulation_initial_name,'_P_',int2str(Dp+Dk),'PE_','Gauss','_N_',int2str(NoiseIn*1e3),'mV (',int2str(q),')');

save(Error_name,'NerrMulti','errMulti');
