% script created by Richard Balson 22/02/2013

% description
% ~~~~~~~~~~~
% this script prints figures to .fig and .pdf if required

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~

% Beginning of script
% ~~~~~~~~~~~~~~~~~~~~~

simulation_name = strcat(Estimation_type,'\',simulation_initial_name,'P',int2str(Dp+Dk),'N',int2str(NoiseIn*1e3),'mV_');

simulation_namez = strcat(simulation_name,'z', int2str(zoomtime),'s_');

% Determine whether figures need to be saved

Image_handling_states_names = {'Vpo' 'Veo' 'Vsio' 'Vfio' 'Zpo' 'Zeo' 'Zsio' 'Zfio' 'Input' 'Excitation' 'Slow_Inhibition' 'Fast_Inhibition'...
    ; 'VpoZ' 'VeoZ' 'VsioZ' 'VfioZ' 'ZpoZ' 'ZeoZ' 'ZsioZ' 'ZfioZ' 'InputZ' 'ExcitationZ' 'Slow_InhibitionZ' 'Fast_InhibitionZ'};

for k = 1:size(Image_handling_states,2) % Determine which states have been plotted and save them as a .fig
    for j = 1:size(Image_handling_states,1)
        
        if (Image_handling_states(j,k) ==1)
            if (j ==1)
                Image_index = NMS(k);
            else
                Image_index = NMSZ(k);
            end
            
            saveas(Image_index,[simulation_name,Image_handling_states_names{j,k},'.fig'],'fig');
            
            if printpdf
                
                print(Image_index, '-dpdf','-painters', '-r2400', [simulation_name,Image_handling_states_names{j,k},'.pdf']);
                
            end
            
        end
        
    end
end

Image_handling_input_names = {'Vp' 'Ve' 'Vsi' 'Vfi'...
    ; 'VpZ' 'VeZ' 'VsiZ' 'VfiZ'};

for k = 1:size(Image_handling_inputs,2) % Determine which states inputs have been ploted and save them as .fig
    for j = 1:size(Image_handling_inputs,1)
        
        if (Image_handling_inputs(j,k) ==1)
            if (j ==1)
                Image_index = NMSI(k);
            else
                Image_index = NMSIZ(k);
            end
            
            saveas(Image_index,[simulation_name,Image_handling_input_names{j,k},'.fig'],'fig');
            
            if printpdf
                
                print(Image_index, '-dpdf','-painters', '-r2400', [simulation_name,Image_handling_input_names{j,k},'.pdf']);
                
            end
            
        end
        
    end
end

Image_handling_multi_names = {'Model_States' 'Model_States_Inputs' 'Model_Parameters' 'Model_Input_and_Parameters';...
    'Model_StatesZ' 'Model_States_InputsZ' 'Model_ParametersZ' 'Model_Input_and_ParametersZ' };

for k = 1:size(Image_handling_multi,2) % Determine which states inputs have been ploted and save them as .fig
    for j = 1:size(Image_handling_multi,1)
        if (((k >2) && Dp+Dk>0) || k <3)
            
            if (Image_handling_multi(j,k) ==1)
                if (j ==1)
                    Image_index = NMM(k);
                else
                    Image_index = NMMZ(k);
                end
                
                saveas(Image_index,[simulation_name,Image_handling_multi_names{j,k},'.fig'],'fig');
                
                if printpdf
                    
                    print(Image_index, '-dpdf','-painters', '-r2400', [simulation_name,Image_handling_multi_names{j,k},'.pdf']);
                    
                end
                
            end
        end
    end
end

if (Estimation_Error && Dp+Dk>0)
    
    Image_index = EstE;
    
    saveas(Image_index,[simulation_name,'Multiple_simulations_Parameters_Error.fig'],'fig');
    
    if printpdf
        
        print(Image_index, '-dpdf','-painters', '-r2400', [simulation_name,'Multiple_simulations_Parameters_Error.pdf']);
        
    end
    
end

