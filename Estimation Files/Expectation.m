% script created by Richard Balson 06/09/2012

% description
% ~~~~~~~~~~~
% This function describes tthe calculation of expected values and covariances using sigma
% points ( Note variances that are described initailly may ned to be added
% to the variances determined using this fucntion

% last edit
% ~~~~~~~~~

% next edit
% ~~~~~~~~~

% Beginning of function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [ExpX ExpY PxxF PxyF PyyF] = Expectation(X, Dx, Y, Dy,kappa)

if nargin == 4
    kappa =0;
end

% Determine element weights
% ~~~~~~~~~~~~~~~~~~~~~~~~~
if kappa > 0
    Weights = [kappa/(Dx+kappa) ones(1,size(X,2)-1)*1/(2*(kappa+Dx))];
    Number_sigma = 2*Dx+1;
else
    Weights = ones(1,size(X,2))*1/(2*(kappa+Dx));
    Number_sigma = 2*Dx;
end

% Determine the expected value for states and output
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

if size(Y,1) == Dy          % Determine how Y matrix orientated
    ExpY = sum(Weights'.*Y')';        % Sum components of each state
else
    ExpY = sum(Weights'.*Y)';
end
if size(X,1) == Dx         % Determine how x matrix orientated
    ExpX= sum(bsxfun(@times,X',Weights'))';
else
    ExpX = sum(bsxfun(@times,X,Weights'))';
end


for i = 1:length(ExpX)
    if size(X,1) == Dx          % Determine how X matrix orientated
        Px(i,:) =  X(i,:) -ExpX(i);
    else 
        Px(i,:) = X(:,i)' - ExpX(i);
    end
end
    if size(Y,1) == Dy          % Determine how X matrix orientated
        Py =  Y -ExpY;
    else 
        Py = Y' - ExpY;
    end

% for i = 1:Number_sigma
%     if size(X,1) == Dx          % Determine how X matrix orientated
%         Px(:,i) =  X(:,i) -ExpX;
%     else 
%         Px(:,i) = X(i,:)' - ExpX;
%     end
%     if size(Y,1) == Dy          % Determine how X matrix orientated
%         Py(:,i) =  Y(:,i) -ExpY;
%     else 
%         Py(:,i) = Y(i,:)' - ExpY;
%     end
% end

PxxF = bsxfun(@times,Px,Weights)*Px';
PxyF = bsxfun(@times,Px,Weights)*Py';
PyyF = (Weights.*Py)*Py';

% for i = 1:Number_sigma
%    if size(X,1) == Dx  
%     PxxT(:,:,i) = Weights(:,i)*(X(:,i)-ExpX)*(X(:,i)-ExpX)';
%     else
%      PxxT(:,:,i) = Weights(:,i)*(X(:,i)'-ExpX)*(X(:,i)'-ExpX)';   
% end
%     if size(Y,1) == Dy          % Determine how X matrix orientated
%         PyyT(:,:,i) =  Weights(:,i).*(Y(:,i) -ExpY)*(Y(:,i) -ExpY)';
%     else 
%         PyyT(:,:,i) = Weights(:,i).*(Y(i,:)' - ExpY)*(Y(:,i)' -ExpY)';
%     end
%     if ((size(X,1)==Dx) && (size(Y,1)==Dy))
%         PxyT(:,:,i) =WeightX(:,i).*(X(:,i) -ExpX)*(Y(:,i) -ExpY)';
%     elseif ((size(X,1)==Dx))
%         PxyT(:,:,i) =WeightX(:,i).*(X(:,i) -ExpX)*(Y(:,i)' -ExpY)';
%     elseif (size(Y,1)==Dy)
%         PxyT(:,:,i) =WeightX(:,i).*(X(:,i)' -ExpX)*(Y(:,i) -ExpY)';
%     else   
%         PxyT(:,:,i) =WeightX(:,i).*(X(:,i)' -ExpX)*(Y(:,i)' -ExpY)';
%     end
% end
% 
% 
% PxxF = sum(PxxT,3);
% PyyF = sum(PyyT,3);
% PxyF = sum(PxyT,3);


% Determine covariance for each state
% ~~~~~~~~~~~~~~~~~~~~~
% XReshape = reshape(X,Dx,1,Number_sigma);
% ExpXReshape = repmat(ExpX,[1 1 Number_sigma]);
% YReshape = reshape(Y,size(Y,1),1,Number_sigma);
% ExpYReshape = repmat(ExpY,[1 1 Number_sigma]);
% WeightsReshape = reshape(Weights,1,1,Number_sigma);
% 
%     PxxF = WeightsReshape*(XReshape-ExpXReshape)*(XReshape-ExpXReshape)';   
%         PyyF =  WeightsReshape.*(YReshape -ExpYReshape)*(YReshape -ExpYReshape)';
%         PxyF =WeightsReshape.*(XReshape -ExpXReshape)*(YReshape -ExpYReshape)';




