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

% Determine whether figures need to be saved

Image_handling_states_names = {'Vpo' 'Veo' 'Vsio' 'Vfio' 'Zpo' 'Zeo' 'Zsio' 'Zfio' 'Input' 'Excitation' 'Slow_Inhibition' 'Fast_Inhibition'...
    ; 'VpoZ' 'VeoZ' 'VsioZ' 'VfioZ' 'ZpoZ' 'ZeoZ' 'ZsioZ' 'ZfioZ' 'InputZ' 'ExcitationZ' 'Slow_InhibitionZ' 'Fast_InhibitionZ'};

for k = Ds+1:size(Image_handling_states,2) % Determine which states have been plotted and save them as a .fig
    
    if (Image_handling_states(1,k))
        
        Image_index = NMSM(k-Ds);
        
        saveas(Image_index,[simulation_name,Image_handling_states_names{1,k},'.fig'],'fig');
        
        if printpdf
            
            print(Image_index, '-dpdf','-painters', '-r2400', [simulation_name,Image_handling_states_names{j,k},'.pdf']);
            
        end
        
    end
    
end

if ((Image_handling_multi(1,3)) && Dp+Dk>0)
    Image_index = NMMM;
    
    saveas(Image_index,[simulation_name,'Multiple_simulations_Parameters.fig'],'fig');
    
    if printpdf
        
        print(Image_index, '-dpdf','-painters', '-r2400', [simulation_name,'Multiple_simulations_Parameters.pdf']);
        
    end
    
end

if (Estimation_Error && Dp+Dk>0)
    
    Image_index = EstEM;
    
    saveas(Image_index,[simulation_name,'Multiple_simulations_Parameters_Error.fig'],'fig');
    
    if printpdf
        
        print(Image_index, '-dpdf','-painters', '-r2400', [simulation_name,'Multiple_simulations_Parameters_Error.pdf']);
        
    end
    
end

