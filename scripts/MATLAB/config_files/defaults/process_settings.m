    % h_data: defines wheter the inference simulation uses ideal values of
    % H or obtains them from a writting opeation. It is defined based on
    % use_write_data
    if use_write_data==1
        h_data='sim';
    else
        h_data='ideal';
    end

    if strcmpi(sim_degradation,'yes')
        if strcmpi(flr_type,'STCKF_ON') || strcmpi(flr_type,'STCKF_OFF') || strcmpi(flr_type,'CORRF_ON') || strcmpi(flr_type,'STCKF_OFF_no_electroformed')
            ratioSTCKF=0.02;  
            deg_folder=strrep(sprintf('%s_%.3f',flr_type,ratioSTCKF),'.','p');
        elseif strcmpi(flr_type,'STCK_MIX')
            deg_folder=strrep(sprintf('%s_%.3f_SA1_%.3f',flr_type,STCKF_total_ratio,ratioSTCKF),'.','p');
        end
    else
    	deg_folder='fresh';
        STCKFratio=[0];
    end
    
    if strcmpi(sim_variability,'yes')
        if strcmpi(variability_mode,'mode1')
            stdImax_array=0;
            stdImin_array=0;          
        elseif strcmpi(variability_mode,'mode2')
            variability_MC=1;
            variabilities=0;      
        end
    else
        stdImax_array=0;
        stdImin_array=0;       
        variability_MC=1;
        variabilities=0;              
    end
    
    if strcmpi(variability_mode,'mode1') && strcmpi(sim_degradation,'no') && strcmpi(timing_stats,'no')
        MC_total=1;   
    end
    
    if strcmpi(dualSide_connect_h,'yes') && strcmpi(dualSide_connect_v,'yes')
        connections='dualSide_connect';
    elseif strcmpi(dualSide_connect_h,'yes') && strcmpi(dualSide_connect_v,'no')
        connections='dualSide_connect_inputs';
    elseif strcmpi(dualSide_connect_h,'no') && strcmpi(dualSide_connect_v,'yes')
    	connections='dualSide_connect_outputs';
    else
    	connections='singleSide_connect';
    end
        
    if strcmpi(correct_rs,'yes')
        if strcmpi(correct_mode,'mode1')
            correct_folder='rs_corrected_mode1';
        elseif strcmpi(correct_mode,'mode2')
            correct_folder='rs_corrected_mode2';
        elseif strcmpi(correct_mode,'mode3')
            correct_folder='rs_corrected_mode3';
        elseif strcmpi(correct_mode,'mode4')
            correct_folder='rs_corrected_mode4';
        elseif strcmpi(correct_mode,'mode5')
            correct_folder='rs_corrected_mode5';
        elseif strcmpi(correct_mode,'mode6')
            correct_folder='rs_corrected_mode6';
        end
    else
        correct_folder='no_correction';
        correct_mode=[];
    end
    
    if size(memdiode_model_array,1)<size(memdiode_model_array,2)
        memdiode_model_array=memdiode_model_array';
    end
    
    [hostname_ret,hostname_name]=system('hostname');
    if strfind(hostname_name,'nanolab5')
        N_Proc=min(max(N_Proc,2),5);
    elseif strfind(hostname_name,'nanolab4')
        N_Proc=min(max(N_Proc,2),5);    
    elseif strfind(hostname_name,'nanolab3')
        N_Proc=min(max(N_Proc,2),15);
    elseif strfind(hostname_name,'nanolab2')
        N_Proc=min(max(N_Proc,2),7);
    elseif strfind(hostname_name,'nanolab1')
        N_Proc=min(max(N_Proc,2),7);
    elseif strfind(hostname_name,'ulabelectron')
        N_Proc=min(max(N_Proc,2),3);
    end
    
    if limit_output==1 && strcmpi(in_polarity,'unipolar')
        polarity=0;
        lims='limited_range';
    elseif limit_output==0  && strcmpi(in_polarity,'unipolar')
        lims='unlimited_range';
        polarity=0;
    elseif limit_output==1  && strcmpi(in_polarity,'bipolar')
        lims='unlimited_range_bipolar';
        polarity=0.5;
    elseif limit_output==0  && strcmpi(in_polarity,'bipolar')
        lims='unlimited_range_bipolar';
        polarity=0.5;
    end
    
    if strcmpi(database,'CIFAR-10')
        database=database;
    elseif strcmpi(database,'SVHN')
        database=database;
    elseif strcmpi(database,'YALE-FACE')
        database=database;    
    elseif strcmpi(database,'YALE-FACE-B')
        database=database;        
    elseif strcmpi(database,'ORL')
        database=database;        
    elseif strcmpi(database,'F-MNIST')
        database=database;
    elseif strcmpi(database,'K-MNIST')
        database=database;    
    else
        database='MNIST';
    end
    
    if strcmpi(learn_algorithm,'trainscg')
        learning_algorithm='Scaled_Conjugate_Gradient';
    else
        learning_algorithm='Scaled_Conjugate_Gradient';
    end
    
    if all_logsig==1
        learning_algorithm=strcat(learning_algorithm,'_allLogsig');
    end
    
    if strcmpi(remap_mode,'mode1')
        remap_folder='remap_mode1';
    elseif strcmpi(remap_mode,'mode2')
        remap_folder='remap_mode2';
    elseif strcmpi(remap_mode,'mode3')
        remap_folder='remap_mode3';
    elseif strcmpi(remap_mode,'mode4')
        remap_folder='remap_mode4';
    elseif strcmpi(remap_mode,'mode5')
        remap_folder='remap_mode5';
    elseif strcmpi(remap_mode,'mode6')
        remap_folder='remap_mode6';
    elseif strcmpi(remap_mode,'mode7')
        remap_folder='remap_mode7';
    elseif strcmpi(remap_mode,'mode8')
        remap_folder='remap_mode8';
    else
        remap_folder='remap_mode1';
    end
    
    if strcmpi(use_ftSize,'yes')
        T=feature_size;
        W=feature_size;
        L=feature_size;
        S=feature_size;
        H=20e-9;
        
        % Capacitance
        % [1] J. Liang, S. Yeh, S. Simon Wong, and H. S. Philip Wong, “Effect of wordline/bitline scaling on the performance, energy consumption, and reliability of cross-point memory array,” ACM J. Emerg. Technol. Comput. Syst., vol. 9, no. 1, pp. 114, 2013.
        e0=8.854e-12;
        epsilon=e0*20;
        epsilon_RRAM=25;
        tox_RRAM=5e-9;
        C_interline=epsilon*L*2.*(0.03.*(W./H)+0.83*(T./H)-0.07*(T./H)*0.222).*(H./S).^1.34;
        C_line2gnd=epsilon*(L/2).*(1.15*(W./H)+2.8*(T./H).^(0.222));
        C_memdiode=feature_size.^2.*epsilon_RRAM*e0./tox_RRAM;

        % Resistance
        % [1] J. Liang, S. Yeh, S. Simon Wong, and H. S. Philip Wong, “Effect of wordline/bitline scaling on the performance, energy consumption, and reliability of cross-point memory array,” ACM J. Emerg. Technol. Comput. Syst., vol. 9, no. 1, pp. 114, 2013.
        p=0.25;
        l_0=39e-9;
        R=0.3;
        d=W-4e-9;
        omega=W-4e-9;
        rho_Cu=1.9e-8;
        Weff=W-4e-9;
        alpha=(l_0./d).*(R./(1-R));

        series_resistance=(L./(Weff.*T)).*rho_Cu.*(((3/4)*(1-p).*(l_0./omega))+3.*(1/3-alpha/2+alpha.^2-alpha.^3.*log(1+(1./alpha))).^(-1));
    end    
    
    total_iterations=max(size(memdiode_model_array))*size(image_size_mat,1)*length(read_voltage)*length(write_voltage)*length(wrFrequency_vec)*length(Dcycle_vec)*length(number_of_images)*length(R_cs)*length(series_resistance)*MC_total*length(STCKFratio)*variability_MC*length(variabilities)*length(stdImin_array)*length(stdImax_array)*length(train_MC_vector)*length(input_vector_freq);
    if total_iterations > 20 && gen_plot==1
        close_after_save=1;
    elseif gen_plot==0
        close_after_save=0;        
    end
    
    if isunix
        sim_folder_base=fullfile(sim_folder_root,'Neural_network_MLP');
    else
        sim_folder_base='D:\NN_RS_UAB\Simulations';
    end
    sim_folder=sim_folder_base;
    folder_index=0;
    while exist(sim_folder,'dir')
        folder_index=folder_index+1;
        sim_folder=strcat(sim_folder_base,'-',num2str(folder_index));
    end
    mkdir(sim_folder);
    
    res_folder_layers='Neural_Network_MLP';  
    
    if strcmpi(neurons,'cmos')
        if strcmpi(senseMode,'TIA')
            sense_folder=sprintf('TIA_gain_%d_Rin_%.2f_cmos',TIA_gain,R_in_TIA);
        elseif strcmpi(senseMode,'Gsense')
            sense_folder=sprintf('Gsense_Rsense_%.2f_cmos',Rsense);
        end
    else
        if strcmpi(senseMode,'TIA')
            sense_folder=sprintf('TIA_gain_%d_Rin_%.2f',TIA_gain,R_in_TIA);
        elseif strcmpi(senseMode,'Gsense')
            sense_folder=sprintf('Gsense_Rsense_%.2f',Rsense);
        end        
    end
    
    noise_folder=sprintf('noise_fmax_%.2e_fmin_%.2e_scale_%.2e',noisefmax,noisefmin,noisescale);

    if strcmpi(train_tool,'MATLAB')
        quantization=missing;
    end
    
DATA_SIM{length(series_resistance),length(read_voltage),max(size(image_size_mat))}=[];