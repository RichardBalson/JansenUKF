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

Sample_inverse = 1/limit; % Used to normalise mean squared error

err = zeros(Dx,1); % Initlaise mean square error
Nerr = zeros(Dp+Dk,1); % Initialise normalised error for initial sum

IndexErr=  MeanCheckTStart*sampling_frequency+1:limit;
IndexErrD = (MeanCheckTStart+EstStart)*sampling_frequency:length(output6);

err(1:Ds) = diag(Sample_inverse*(z(IndexErrD,:)' - X(1:Ds,IndexErr))*(z(IndexErrD,:)' - X(1:Ds,IndexErr))');

% Determine mean square error for estimated model output
% ~~~~~~~~~~~~~~~

Oerr = Sample_inverse*(Y(IndexErr)' - ExpY(1,IndexErr))*(Y(IndexErr)' - ExpY(1,IndexErr))';

% Determine mean square error for parameter estimates

if (Dp)>0
    err(Ds+Dk+1:Ds+Dk+Dp) = diag(Sample_inverse*(MVI(IndexErrD,1:Dp)' - X(Ds+Dk+1:Ds+Dk+Dp,IndexErr))*(MVI(IndexErrD,1:Dp)' - X(Ds+Dk+1:Ds+Dk+Dp,IndexErr))');
end

if Dk >0
    err(Ds+1) = Sample_inverse*(ones(1,length(IndexErrD))*Input_mean-X(Ds+1,IndexErr))*(ones(1,length(IndexErrD))*Input_mean - X(Ds+1,IndexErr))';
end

% Determine normalised error

if (Dp)>0
    New = MVI;
    New(New==0) = eps;
    Nerr(1:Dp) = Sample_inverse*diag((MVI(IndexErrD,1:Dp)' - X(Ds+Dk+1:Ds+Dk+Dp,IndexErr))./New(IndexErrD,1:Dp)'*(MVI(IndexErrD,1:Dp)' - X(Ds+Dk+1:Ds+Dk+Dp,IndexErr))');
        PercErrEstimate = abs(X(Ds+Dk+1:Ds+Dk+Dp,1:end-1)-MVI(EstStart*sampling_frequency:end,1:Dp)')./MVI(EstStart*sampling_frequency:end,1:Dp)'*100;
end

if Dk >0
    Nerr(Dp+1) = Sample_inverse*diag((Input_mean-X(Ds+1,IndexErr))/Input_mean*(Input_mean-X(Ds+1,IndexErr))');
        PercErrEstimate(Dp+1,:) = abs(X(Ds+Dk,1:end-1)-Input_mean)/Input_mean*100;
end


