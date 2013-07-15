 % script created by Richard Balson 14/10/2012

% description
% ~~~~~~~~~~~
% This function describes the Wendling neural mass model with 8 differential equations and its
% observation function, the simplified version

% last edit
% ~~~~~~~~~

% Yout function defined correctly

% next edit
% ~~~~~~~~~

% Beginning of function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [Xout Yout] = JNM(Xin,dt,CoV, gain, tcon,C)
Xout = Xin;
zs1 =  Xin(2,:)-Xin(3,:);
[Xout(1,:) Xout(4,:)] = PSPkernel([Xin(1,:); Xin(4,:)],dt,gain(1,:),tcon(1),sigmoid(zs1));
zs2 = C(1)*Xin(1,:);
[Xout(2,:) Xout(5,:)] = PSPkernel([Xin(2,:); Xin(5,:)],dt,gain(1,:),tcon(1),CoV +C(2)*sigmoid(zs2));
zs3 = C(3)*Xin(1,:);
[Xout(3,:) Xout(6,:)] = PSPkernel([Xin(3,:); Xin(6,:)],dt,gain(2,:),tcon(2),C(4)*sigmoid(zs3));

Yout = Xout(2,:)-Xout(3,:); % Observation function

end

