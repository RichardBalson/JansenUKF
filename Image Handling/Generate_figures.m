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

% Plot of output
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

maxlimit = round(sampling_frequency/10):length(check);

% % Plot of state 1
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Error_multiplier = 1;

t = linspace(0,dt*length(check),length(check));
tz = linspace(tstart,tstart+zoomtime,zoomtime*sampling_frequency+1);

index = {1:length(check) tstart*sampling_frequency+1:(tstart+zoomtime)*sampling_frequency+1 ...
    EstStart_Sample:length(output6) round((EstStart+tstart)*sampling_frequency)+1:round((EstStart+tstart+zoomtime)*sampling_frequency)+1};

if (Image_handling_model_output(1,1))
    
    legendT = {'Noisy Obs.','Estimated Obs.','Obs.'};
    NMO=state_figure('Jansen Neural Mass Output','Obs',fig_settings,t,[Y(index{1})'; X(2,index{1})-X(3,index{1}); check(index{1})'],legendT,[]);
    
end

if (Image_handling_model_output(2,1))
    
    legendT = {'Noisy Obs.','Estimated Obs.','Obs.'};
    NMO=state_figure('Jansen Neural Mass Output zoomtimeed In','Obs',fig_settings,tz,[Y(index{2})'; X(2,index{2})-X(3,index{2}); check(index{2})'],legendT,[]);
    
end
legendT = {'Sim. v_{p0}','Est. v_{p0}','Std. Dev.';...
    'Sim. v_{p1}','Est. v_{p1}','Std. Dev.';...
    'Sim. v_{p2}','Est. v_{p2}','Std. Dev.';...
    'Sim. z_{p0}','Est. z_{p0}','Std. Dev.';...
    'Sim. z_{p1}','Est. z_{p1}','Std. Dev.';...
    'Sim. z_{p2}','Est. z_{p2}','Std. Dev.'};

fig_name = {'Neural Mass State1 (v_{p0})';...
    'Neural Mass State2 (v_{p1})';...
    'Neural Mass State3 (v_{p2})';...
    'Neural Mass State4 (z_{p0})';...
    'Neural Mass State5 (z_{p1})';...
    'Neural Mass State6 (z_{p2})'};


for k = 1:Ds
    
    if (Image_handling_states(1,k))
        % Plot of state Not zoomtimeed in
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (plot_uncertainty)
        erfn = X(k,index{1})-squeeze(sqrt(Pxx(k,k,index{1})*Error_multiplier))';
        erfp = X(k,index{1})+squeeze(sqrt(Pxx(k,k,index{1})*Error_multiplier))';
        else
            erfn =[]; erfp =[];
        end
        NMS(k)=state_figure(fig_name{k},'State',fig_settings,t,[z(index{3},k)';X(k,index{1})],legendT(k,:),[erfn;erfp]);
        
    end
end
clear erfn erfp
for k = 1:Ds
    
    if (Image_handling_states(2,k))
        % Plot of state Not zoomtimeed in
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (plot_uncertainty)
        erfn = X(k,index{2})-squeeze(sqrt(Pxx(k,k,index{2})*Error_multiplier))';
        erfp = X(k,index{2})+squeeze(sqrt(Pxx(k,k,index{2})*Error_multiplier))';
        else
            erfn =[]; erfp =[];
        end
        NMSZ(k)=state_figure(fig_name{k},'State',fig_settings,tz,[z(index{4},k)';X(k,index{2})],legendT(k,:),[erfn;erfp]);
        
    end
end
clear erfn erfp
% Plot input state
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~

if Dk
    if (Image_handling_states(1,Ds+1))
        if (plot_uncertainty) 
        erfn = X(Ds+Dk,index{1})-squeeze(sqrt(Pxx(Ds+Dk,Ds+Dk,index{1})*Error_multiplier))';
        erfp = X(Ds+Dk,index{1})+squeeze(sqrt(Pxx(Ds+Dk,Ds+Dk,index{1})*Error_multiplier))';
        else
            erfn =[]; erfp =[];
        end
        NMS(Ds+Dk)=state_figureInput('Neural Mass (Input)','State',fig_settings,t,[meanf*ones(1,length(t));X(Ds+Dk,index{1})],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        
    end
    
    if (Image_handling_states(2,Ds+1))
      if (plot_uncertainty)  
        erfn = X(Ds+Dk,index{2})-squeeze(sqrt(Pxx(Ds+Dk,Ds+Dk,index{2})*Error_multiplier))';
        erfp = X(Ds+Dk,index{2})+squeeze(sqrt(Pxx(Ds+Dk,Ds+Dk,index{2})*Error_multiplier))';
        else
            erfn =[]; erfp =[];
        end
        NMSZ(Ds+Dk)=state_figureInput('Neural Mass (Input)','State',fig_settings,tz,[meanf*ones(1,length(tz));X(Ds+Dk,index{2})],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        
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
        if (Image_handling_states(1,Ds+1+k))
           if (plot_uncertainty)  
            erfn = X(Ds+Dk+k,index{1})-squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))';
            erfp = X(Ds+Dk+k,index{1})+squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))';
            else
            erfn =[]; erfp =[];
        end
            NMS(Ds+Dk+k)=state_figure(fig_name{k},'State',fig_settings,t,[MVI(index{3},k)';X(Ds+Dk+k,index{1})],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        end
        if (Image_handling_states(2,Ds+1+k))
            if (plot_uncertainty) 
            erfn = X(Ds+Dk+k,index{2})-squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{2})*Error_multiplier))';
            erfp = X(Ds+Dk+k,index{2})+squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{2})*Error_multiplier))';
            else
            erfn =[]; erfp =[];
        end
            NMSZ(Ds+Dk+k)=state_figure(fig_name{k},'State',fig_settings,tz,[MVI(index{4},k)';X(Ds+Dk+k,index{2})],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        end
    end
end
clear erfn erfp

%%

% Plot inputs to each neural population
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig_name = {'Neural Mass Membrane Potential (v_{p})';...
    'Neural Mass Membrane Potential (v_{e})';...
    'Neural Mass Membrane Potential (v_{s})'};
if any(Image_handling_inputs(1,:))
    Input(:,:,1) = [Y(index{1})'; ExpY(index{1})];
    Input(:,:,2) = [C(1)*z(index{3},1)'; C(1)*X(1,index{1})];
    Input(:,:,3) = [C(3)*z(index{3},1)'; C(3)*X(1,index{1})];
end
for k = Ds/2
    if (Image_handling_inputs(1,k))
        NMSI(k)=state_figure(fig_name{k},'State',fig_settings,t,Input(:,:,k),{'Sim. Input','Est. Input'},[]);
    end
end
if any(Image_handling_inputs(2,:))
    clear Input
    Input(:,:,1) = [Y(index{2})'; ExpY(index{2})];
    Input(:,:,2) = [C(1)*z(index{4},1)'; C(1)*X(1,index{2})];
    Input(:,:,3) = [C(3)*z(index{4},1)'; C(3)*X(1,index{2})];
end
for k = Ds/2
    if (Image_handling_inputs(2,k))
        NMSIZ(k)=state_figure(fig_name{k},'State',fig_settings,tz,Input(:,:,k),{'Sim. Input','Est. Input'},[]);
    end
end
clear Input

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all aggregate membrane potentials on a single plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~

if (Image_handling_multi(1,1) ==1) % Plot all model states
    RowP = Ds/2; % Set number of rows for multi plot
    ColP =2; % Set number of columns for multi plot
    fig_name = 'All Fast States';
    PlotsPerFig =2;
    for k = 1:Ds
        input((k-1)*PlotsPerFig+1:k*PlotsPerFig,:) =[z(index{3},k)';X(k,index{1})];
        if plot_uncertainty
          erfn(k,:) = X(k,index{1})-squeeze(sqrt(Pxx(k,k,index{1})*Error_multiplier))';
          erfp(k,:) = X(k,index{1})+squeeze(sqrt(Pxx(k,k,index{1})*Error_multiplier))';
        else erfn =[]; erfp =[];
        end
    end
    namesS = {'v_{po} (mV)','v_{eo} (mV)','v_{sio} (mV)','Z_{po} (mV)','Z_{eo} (mV)','Z_{sio} (mV)'};
    NMM(1) = state_figure_multi(fig_name,'State',fig_settings,t,input,{'Sim.','Est.','Std. Dev'},[erfn;erfp],RowP,ColP,namesS,PlotsPerFig);
    clear input erfn erfp namesS
end
    
if (Image_handling_multi(2,1) ==1) % Plot all model states
       RowP = Ds/2; % Set number of rows for multi plot
    ColP =2; % Set number of columns for multi plot
    fig_name = 'All Fast States';
    PlotsPerFig =2;
    for k = 1:Ds
        input((k-1)*PlotsPerFig+1:k*PlotsPerFig,:) =[z(index{4},k)';X(k,index{2})];
        if plot_uncertainty
          erfn(k,:) = X(k,index{2})-squeeze(sqrt(Pxx(k,k,index{2})*Error_multiplier))';
          erfp(k,:) = X(k,index{2})+squeeze(sqrt(Pxx(k,k,index{2})*Error_multiplier))';
        else erfn =[]; erfp =[];
        end
    end
    namesS = {'v_{po} (mV)','v_{eo} (mV)','v_{sio} (mV)','Z_{po} (mV)','Z_{eo} (mV)','Z_{sio} (mV)'};
    NMMZ(1) = state_figure_multi(fig_name,'State',fig_settings,tz,input,{'Sim.','Est.','Std. Dev'},[erfn;erfp],RowP,ColP,namesS,PlotsPerFig);
    clear input erfn erfp namesS
end

if (Image_handling_multi(1,2) ==1) % Plot all state inputs
    % Plot all aggregate membrane potentials on a single plot
    RowP = 2;
    ColP = 2;
    PlotsPerFig =2;
    input = [Y(index{1})';ExpY(index{1});...
            C(1)*z(index{3},1)';C(1)*X(1,index{1});...
            C(3)*z(index{3},1)';C(3)*X(1,index{1})];
        fig_name = 'Neural Mass State Inputs';
    namesS = {'v_{p} (mV)','v_{e} (mV)','v_{si} (mV)'};
    NMM(2) = state_figure_multi(fig_name,'State',fig_settings,t,input,{'Sim.','Est.'},[],RowP,ColP,namesS,PlotsPerFig);
    clear input erfn erfp namesS
end

if (Image_handling_multi(2,2) ==1) % % Plot all aggregate membrane potentials on a single plot
    RowP = 2;
    ColP = 2;
    PlotsPerFig =2;
    input = [Y(index{2})';ExpY(index{2});...
            C(1)*z(index{4},1)';C(1)*X(1,index{2});...
            C(3)*z(index{4},1)';C(3)*X(1,index{2});...
            (C(5)*z(index{4},1)-C(6)*z(index{4},3))';C(5)*X(1,index{2})-C(6)*X(3,index{2})];
        fig_name = 'Neural Mass State Inputs';
    namesS = {'v_{p} (mV)','v_{e} (mV)','v_{si} (mV)'};
    NMMZ(2) = state_figure_multi(fig_name,'State',fig_settings,tz,input,{'Sim.','Est.'},[],RowP,ColP,namesS,PlotsPerFig);
    clear input erfn erfp namesS    
end

if ((Image_handling_multi(1,3) ==1) &&(Dp+Dk >0)) % Plot model parameters estimated
    RowP =1;
    ColP=1;
    if (Dp+Dk>2)
        RowP =2;
    end
    if (Dp+Dk>1)
        ColP =2;
    end
    PlotsPerFig =2;
    for k =1:Dp
    input((k-1)*PlotsPerFig+1:k*PlotsPerFig,:) = [MVI(index{3},k)';X(Ds+Dk+k,index{1})];
        if plot_uncertainty
          erfn(k,:) = X(Ds+Dk+k,index{1})-squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))';
          erfp(k,:) = X(Ds+Dk+k,index{1})+squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))';
        else erfn =[]; erfp =[];
        end
    end
    if Dk 
        input(Dp*PlotsPerFig+1:(Dp+1)*PlotsPerFig,:)= [meanf*ones(length(t),1)';X(Ds+1,index{1})];
        if plot_uncertainty
          erfn(Dp+1,:) = X(Ds+1,index{1})-squeeze(sqrt(Pxx(Ds+1,Ds+1,index{1})*Error_multiplier))';
          erfp(Dp+1,:) = X(Ds+1,index{1})+squeeze(sqrt(Pxx(Ds+1,Ds+1,index{1})*Error_multiplier))';
          else erfn =[]; erfp =[];
        end
    end
        fig_name = 'Neural Mass Parameters';
    namesS = {'G_{p} (mV)','G_{s} (mV)','Input (Hz)'};
    NMM(3) = state_figure_multi(fig_name,'State',fig_settings,t,input,{'Sim.','Est.'},[erfn;erfp],RowP,ColP,namesS,PlotsPerFig);
    clear input erfn erfp namesS    
    
end

if ((Image_handling_multi(2,3) ==1) && Dp+Dk>0) % Plot model parameters estimated
    RowP =1;
    ColP=1;
    if (Dp+Dk>2)
        RowP =2;
    end
    if (Dp+Dk>1)
        ColP =2;
    end
    PlotsPerFig =2;
    for k =1:Dp
    input((k-1)*PlotsPerFig+1:k*PlotsPerFig,:) = [MVI(index{4},k)';X(Ds+Dk+k,index{2})];
        if plot_uncertainty
          erfn(k,:) = X(Ds+Dk+k,index{2})-squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{2})*Error_multiplier))';
          erfp(k,:) = X(Ds+Dk+k,index{2})+squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{2})*Error_multiplier))';
        else erfn =[]; erfp =[];
        end
    end
    if Dk 
        input(Dp*PlotsPerFig+1:(Dp+1)*PlotsPerFig,:)= [meanf*ones(length(tz),1)';X(Ds+1,index{2})];
        if plot_uncertainty
          erfn(Dp+1,:) = X(Ds+1,index{2})-squeeze(sqrt(Pxx(Ds+1,Ds+1,index{2})*Error_multiplier))';
          erfp(Dp+1,:) = X(Ds+1,index{2})+squeeze(sqrt(Pxx(Ds+1,Ds+1,index{2})*Error_multiplier))';
          else erfn =[]; erfp =[];
        end
    end
        fig_name = 'Neural Mass Parameters';
    namesS = {'G_{p} (mV)','G_{s} (mV)','Input (Hz)'};
    NMM(3) = state_figure_multi(fig_name,'State',fig_settings,tz,input,{'Sim.','Est.'},[erfn;erfp],RowP,ColP,namesS,PlotsPerFig);
    clear input erfn erfp namesS  
end


%%

if (Image_handling_firing_rates(1) ==1)
    % Plot of state 1
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RowP =3;
    ColP =1;
    PlotsPerFig =2;
    fig_name = 'Neural Mass Pyramidal Neuron Input, Output and Firing Rate';
    input = [Y(index{1})';ExpY(index{1});...
            Sigmoid(Y(index{1}))';Sigmoid(ExpY(index{1}));...
            z(index{3},1)'; X(1,index{1})];
    namesS = {'Pyramidal Input (mV)' 'Pyramidal Firing Rate (Hz)' 'Pyramidal Output (mV)'};
        NMFR(1) = state_figure_multi(fig_name,'State',fig_settings,t,input,{'Sim.','Est.'},[],RowP,ColP,namesS,PlotsPerFig);
    clear input namesS
end

if (Image_handling_firing_rates(2) ==1)
    
    % Plot of state 2
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RowP =3;
       ColP =1;
    PlotsPerFig =2;
    fig_name = 'Neural Mass Excitatory Neuron Input, Output and Firing Rate';
    input = [C(1)*z(index{3},1)';C(1)*X(1,index{1});...
            Sigmoid(C(1)*z(index{3},1))';Sigmoid(C(1)*X(1,index{1}));...
            z(index{3},2)';X(2,index{1})];
    namesS = {'Excitatory Input (mV)' 'Excitatory Firing Rate (Hz)' 'Excitatory Output (mV)'};
    NMFR(2) = state_figure_multi(fig_name,'State',fig_settings,t,input,{'Sim.','Est.'},[],RowP,ColP,namesS,PlotsPerFig);
    clear input namesS
end

if (Image_handling_firing_rates(3) ==1)
    
    % Plot of state 3
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RowP =3;
    ColP=1;
    PlotsPerFig =2;
    fig_name = 'Neural Mass Excitatory Neuron Input, Output and Firing Rate';
    input = [C(3)*z(index{3},1)';C(3)*X(1,index{1});...
            Sigmoid(C(3)*z(index{3},1))';Sigmoid(C(3)*X(1,index{1}));...
            z(index{3},3)';X(3,index{1})];
    fig_name = 'Neural Mass Slow Inh. Neuron Input, Output and Firing Rate';
%     MassP = {C(3)*z(index{3},1); Sigmoid(C(3)*z(index{3},1)); z(index{3},3)};
%     MassPE = {C(3)*X(1,index{3}); Sigmoid(C(3)*X(1,index{3})); X(3,index{3})};
    namesS = {'Slow Inh. Input (mV)' 'Slow Inh. Firing Rate (Hz)' 'Slow Inh. Output (mV)'};
        NMFR(3) = state_figure_multi(fig_name,'State',fig_settings,t,input,{'Sim.','Est.'},[],RowP,ColP,namesS,PlotsPerFig);
    clear input namesS
end

if (Estimation_Error && Dp+Dk>0)
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
    input(k,:) = PercErrEstimate(k,:);
        if plot_uncertainty
          erfn(k,:) = PercErrEstimate(k,:)-squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))'./MVI(index{3},k)'*100;
          erfp(k,:) = PercErrEstimate(k,:)+squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))'./MVI(index{3},k)'*100;
        else erfn =[]; erfp =[];
        end
    end
    if Dk 
        input(Dp+1,:) = PercErrEstimate(Dp+1,:);
        if plot_uncertainty
          erfn(Dp+1,:) = PercErrEstimate(Dp+1,:)-squeeze(sqrt(Pxx(Ds+1,Ds+1,index{1})*Error_multiplier))'./Input_mean*100;
          erfp(Dp+1,:) = PercErrEstimate(Dp+1,:)+squeeze(sqrt(Pxx(Ds+1,Ds+1,index{1})*Error_multiplier))'./Input_mean*100;
          else erfn =[]; erfp =[];
        end
    end
        fig_name = 'Estimation Error';
    namesS = {'G_{p} (% Error)','G_{s} (% Error)','Input (% Error)'};
    EstE(1) = state_figure_multi(fig_name,'State',fig_settings,t,input,{'Estimation Error'},[erfn;erfp],RowP,ColP,namesS,PlotsPerFig);
    clear input erfn erfp namesS    
%     EstE=state_figure(fig_name,'Obs',fig_settings,t,PercErrEstimate,legen
%     dT,[]);
end

