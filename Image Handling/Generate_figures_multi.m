% script created by Richard Balson 20/01/2013

% description
% ~~~~~~~~~~~
% this script generates figures for UKF estimation results for the Document
% for the Boon group

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~

% At line 184 editing multi figure plots

% Beginning of script
% ~~~~~~~~~~~~~~~~~~~~~

% Generic Figure handling parameters

fig_settings.label_fontsize = 10;            % point
fig_settings.tick_fontsize = 8;              % point
fig_settings.legend_fontsize = 10;
fig_settings.left_pos = 5;               % cms
fig_settings.bottom_pos = 5;             % cms
fig_settings.font_type = 'Arial';
fig_settings.dirname = 'Results';              % default directory for figure files
fig_settings.scale =0.5;

% % Plot of state 1
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Error_multiplier = 1;
scale = 0.5;

t = linspace(0,dt*Decimate*size(X_Multi,1),size(X_Multi,1));

index = {1:size(X_Multi,1) (EstStart)*sampling_frequency+1:Decimate:length(output6)};

% Plot input state
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~

if Dk
    if (Image_handling_states(1,9))
        if (plot_uncertainty)
            erfn = squeeze(X_Multi(:,1,:)-sqrt(Pxx_Multi(:,1,index{1})*Error_multiplier))';
            erfp = squeeze(X_Multi(:,1,:)+sqrt(Pxx_Multi(:,1,index{1})*Error_multiplier))';
        else
            erfn =[]; erfp =[];
        end
        NMSM(Dk)=state_figureInput('Neural Mass (Input)','State',fig_settings,t,[meanf*ones(length(t),1);squeeze(X_Multi(index{1},Ds+Dk,:))],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        
    end
end
clear erfn erfp
% Determine whether to plot parameter estimates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~

fig_name = {'Neural Mass Exc. Gain (G_{p})';...
    'Neural Mass Slow Inh. Gain (G_{s})'};
if Dp>0
    for k = 1:Dp
        if (Image_handling_states(1,Ds+k+1))
            if (plot_uncertainty)
                erfn = squeeze(X_Multi(:,k+Dk,:)-sqrt(Pxx_Multi(:,k+Dk,:)*Error_multiplier))';
                erfp = squeeze(X_Multi(:,k+Dk,:)+sqrt(Pxx_Multi(:,k+Dk,:)*Error_multiplier))';
            else
                erfn =[]; erfp =[];
            end
            NMSM(Dk+k-1)=state_figure(fig_name{k},'State',fig_settings,t,[MVI(index{2},k)';squeeze(X_Multi(:,Dk+k,:))],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        end
    end
end
clear erfn erfp input

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all aggregate membrane potentials on a single plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~

if ((Image_handling_multi(1,3) ==1) &&(Dp+Dk >0)) % Plot model parameters estimated
    RowP =1;
    ColP=1;
    if (Dp+Dk>2)
        RowP =2;
    end
    if (Dp+Dk>1)
        ColP =2;
    end
    PlotsPerFig =1 + Simulation_number;
    for k =1:Dp
        input(:,:,k) = [MVI(index{2},k)';squeeze(X_Multi(:,Dk+k,:))'];
        if plot_uncertaintyMulti
                erfn(:,:,k) = squeeze(X_Multi(:,k+Dk,:)-sqrt(Pxx_Multi(:,k+Dk,:)*Error_multiplier))';
                erfp(:,:,k) = squeeze(X_Multi(:,k+Dk,:)+sqrt(Pxx_Multi(:,k+Dk,:)*Error_multiplier))';
        else erfn =[]; erfp =[];
        end
    end
    if Dk
        input(:,:,Dp+1)= [meanf*ones(length(t),1)';squeeze(X_Multi(:,1,:))'];
        if plot_uncertaintyMulti
            erfn(:,:,Dp+1) = X_Multi(:,1,:)-squeeze(sqrt(Pxx(1,1,:)*Error_multiplier))';
            erfp(:,:,Dp+1) = X_Multi(:,1,:)+squeeze(sqrt(Pxx(1,1,:)*Error_multiplier))';
        else erfn =[]; erfp =[];
        end
    end
    fig_name = 'Neural Mass Parameters';
    namesS = {'G_{p} (mV)','G_{s} (mV)','Input (Hz)'};
    NMMM = state_figure_multi_sim(fig_name,'Multi',fig_settings,t,input,{'Sim.','Est.'},erfn,erfp,RowP,ColP,namesS,PlotsPerFig);
    clear input erfn erfp namesS
    
end

if Estimation_Error
    RowP =1;
    ColP=1;
    if (Dp+Dk>2)
        RowP =2;
    end
    if (Dp+Dk>1)
        ColP =2;
    end
    PlotsPerFig =1;
    for k =1:Dp
        input(:,:,k)= squeeze(PercErrEstimateMulti(k+Dk,:,:))';
        if plot_uncertaintyMulti
            erfn(:,:,k) = bsxfun(@rdivide,squeeze(PercErrEstimateMulti(k+Dk,:,:))'-squeeze(sqrt(Pxx_Multi(:,k+Dk,:))*Error_multiplier)',MVI(index{2},k)'*100);
            erfp(:,:,k) = bsxfun(@rdivide,squeeze(PercErrEstimateMulti(k+Dk,:,:))'+squeeze(sqrt(Pxx_Multi(:,k+Dk,:))*Error_multiplier)',MVI(index{2},k)'*100);
        else erfn =[]; erfp =[];
        end
    end
    if Dk
        input(:,:,Dp+1) = squeeze(PercErrEstimateMulti(1,:,:))';
        if plot_uncertaintyMulti
            erfn(:,:,Dp+1) = squeeze(PercErrEstimateMulti(1,:,:))'-squeeze(sqrt(Pxx_Multi(:,1,:))*Error_multiplier)'./Input_mean*100;
            erfp(:,:,Dp+1) = squeeze(PercErrEstimateMulti(1,:,:))'+squeeze(sqrt(Pxx_Multi(:,1,:))*Error_multiplier)'./Input_mean*100;
        else erfn =[]; erfp =[];
        end
    end
    fig_name = 'Estimation Error';
    namesS = {'G_{p} (% Error)','G_{s} (% Error)','Input (% Error)'};
    EstE = state_figure_multi_sim(fig_name,'State',fig_settings,t,input,{'Estimation Error'},erfn,erfp,RowP,ColP,namesS,PlotsPerFig);
    clear input erfn erfp namesS
    %     EstE=state_figure(fig_name,'Obs',fig_settings,t,PercErrEstimate,legen
    %     dT,[]);
end


