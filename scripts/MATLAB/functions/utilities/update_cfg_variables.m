function [cfg_variables_new] = update_cfg_variables(cfg_variables,results_folder,varargin)

    mode='complete';
    if mod(nargin,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-2
            if strcmp(varargin{i},'mode')
                mode=varargin{i+1};
            end
        end
    end
    
    cfg_vars_to_remove={'MC_total','wv_dir','DATA_SIM','debug_mode','folder_index',...
                        'varargin','test_write','hostname_name','hostname_ret','hspice_dir',...
                        'LICENSE_address','N_Proc','pos_neg','finesim_dir','print_I_tolerance',...
                        'print_V_tolerance','remove_RAW_after_sim','sim_folder','sim_folder_base',...
                        'sim_folder_root','project_directory','gen_plot','run_sim','close_after_save',...
                        'save_all','save_h','save_imemd','save_inner_neuron','save_vpn','timeSTEP_sim_memdiode',...
                        'save_write_signals','deg_folder','remap_folder','res_folder_layers','sense_folder',...
                        'simulation_read_time','simulation_write_time','simulation_WR_margin_time','noise_folder',...
                        'total_iterations','ratioSTCKF','proc_settings','memdiode_model'};
    
    for cfg_vars_to_remove_i=1:length(cfg_vars_to_remove)
        if isfield(cfg_variables,cfg_vars_to_remove{cfg_vars_to_remove_i})
            cfg_variables = rmfield(cfg_variables,cfg_vars_to_remove{cfg_vars_to_remove_i});   
        end
    end
    
    fields_cfg_variables=fields(cfg_variables);
    Index = find(contains(fields_cfg_variables,'hidden_layer_'));
    for index_i=1:length(Index)
        value=sscanf(fields_cfg_variables{Index(index_i),1},'hidden_layer_%d');
        if value>6
            cfg_variables=rmfield(cfg_variables,fields_cfg_variables{Index(index_i)},1);
        end
    end

    if strcmpi(mode,'complete')
        result_info=split(results_folder,'/');
        for j=1:length(result_info)
            str=result_info{j,1};

            % Memdiode model
            if contains(str,'_Hmin_') && contains(str,'-Hmax_')
                aux_str=split(str,'_Hmin_');
                cfg_variables.memdiode_model_array=aux_str{1,1};
            end

            % Series resistance and R_cs
            if contains(str,'Rs') && contains(str,'Rcs')
                values=sscanf(strrep(str,'p','.'),'Rs_%f_Rcs_%f');
                cfg_variables.series_resistance=values(1,1);
                cfg_variables.R_cs=values(2,1);
            end

            % Write frequency, rise time of the write signal, Duty cycle of the
            % programming pulse and duty cycle of the verify pulse
            if contains(str,'freq') && contains(str,'tr') && contains(str,'DC') && contains(str,'DCv')
                values=sscanf(str,'freq_%f_tr_%f_DC_%f_DCv_%f');
                cfg_variables.wrFrequency_vec=values(1,1);
                cfg_variables.trise=values(2,1);
                cfg_variables.Dcycle_vec=values(3,1);
                cfg_variables.Dcycle_v=values(4,1);
            end

            % Vread
            if contains(str,'Vread_')
                values=sscanf(strrep(str,'p','.'),'Vread_%f');
                cfg_variables.read_voltage=values(1,1);
            end

            % Vwrite
            if contains(str,'Vwrite')
                values=sscanf(strrep(str,'p','.'),'Vwrite_%f');
                cfg_variables.write_voltage=values(1,1);
            end

            % Numer of images
            if contains(str,'N_')
                if contains(str(1:2),'N_')
                    values=sscanf(str,'N_%d');
                    cfg_variables.number_of_images=values(1,1);
                end
            end

            % Image resolution
            if contains(str,'_px') 
                values=sscanf(str,'%dx%d_px');
                cfg_variables.image_size_mat={[values(1,1) values(2,1)]};
                cfg_variables.neu_per_layer={[values(1,1)*values(2,1) 10]};
            end

            % Particiones
            if contains(str,'Partitions') 
                num_dim=length(strfind(str,'_'));
                for num_dim_i=1:num_dim
                    if num_dim_i==1
                        reference_string='Partitions_%d_';
                    elseif num_dim_i<num_dim
                        reference_string=strcat(reference_string,'%d_');
                    else
                        reference_string=strcat(reference_string,'%d');
                    end
                end
                values=sscanf(str,reference_string);
                cfg_variables.partitions_N={values'};
            end

            if contains(str,'STCKF_')
                if strcmpi(cfg_variables.sim_degradation,'yes')
                    if strcmpi(cfg_variables.flr_type,'STCKF_ON')
                        values=sscanf(strrep(str,'p','.'),'STCKF_ON_%f');
                        cfg_variables.STCKF_total_ratio=missing;
                    elseif strcmpi(cfg_variables.flr_type,'STCKF_OFF')
                        values=sscanf(strrep(str,'p','.'),'STCKF_OFF_%f');
                        cfg_variables.STCKF_total_ratio=missing;
                    elseif strcmpi(cfg_variables.flr_type,'STCKF_OFF_no_electroformed')
                        values=sscanf(strrep(str,'p','.'),'STCKF_OFF_no_electroformed_%f');
                        cfg_variables.STCKF_total_ratio=missing;
                    elseif strcmpi(cfg_variables.flr_type,'CORRF_ON')
                        values=sscanf(strrep(str,'p','.'),'CORRF_ON_%f');
                        cfg_variables.STCKF_total_ratio=missing;
                    elseif strcmpi(cfg_variables.flr_type,'STCKF_MIX')
                        values=sscanf(strrep(str,'.','p'),'STCKF_MIX_%f_SA1_%f');
                        cfg_variables.STCKFratio=values(1,1);                    
                    end
                    cfg_variables.STCKFratio=values(1,1);
                else
                    cfg_variables.STCKFratio=missing;
                    cfg_variables.STCKF_total_ratio=missing;
                    cfg_variables.flr_type=missing;
                    cfg_variables.reliability_MC=missing;
                end    
            end

            if contains(str,'variability_')
                values=sscanf(strrep(str,'p','.'),'fresh_variability_%f');
                cfg_variables.sim_variability='yes';
                cfg_variables.variabilities=values;
                cfg_variables.stdImin_array=missing;
                cfg_variables.stdImax_array=missing; 
                cfg_variables.variability_mode='mode1';
            elseif strcmpi(cfg_variables.variability_mode,'mode2')
                values=sscanf(strrep(str,'p','.'),'stdImin_%f_stdImax_%f');
                cfg_variables.sim_variability='yes';
                cfg_variables.variabilities=missing;
                cfg_variables.stdImin_array=values(1,1);
                cfg_variables.stdImax_array=values(2,1);   
                cfg_variables.variability_mode='mode2';                     
            end

            if ~contains(results_folder,'variability_') && ~contains(results_folder,'stdImin_') && ~contains(results_folder,'stdImax_')
                cfg_variables.variability_mode=[];
                cfg_variables.variabilities=missing;
                cfg_variables.stdImin_array=missing;
                cfg_variables.stdImax_array=missing;
                cfg_variables.sim_variability='no';
            end

            if ~contains(results_folder,'STCKF_') && ~contains(results_folder,'CORRF_') 
                cfg_variables.STCKFratio=missing;
                cfg_variables.STCKF_total_ratio=missing;
                cfg_variables.flr_type=[];
                cfg_variables.reliability_MC=missing;
                cfg_variables.sim_degradation='no';
            end

            if contains(str,'MLP_')
                if ~isfield(cfg_variables,'neu_per_layer')
                    cfg_variables.neu_per_layer=sscanf(str,'MLP_%dby%d_closed_loop');
                end
            end 

            if contains(str,'STCKF_') && contains(str,'CORRF_') 
                cfg_variables.reliability_MC=sscanf(result_info{j+4,1},'MC= %d');
            else
                cfg_variables.reliability_MC=missing;
            end
        end

        if exist(fullfile(cfg_variables.results_folder, 'calc_pwr.mat'),'file')
            pwr_metrics=load(fullfile(cfg_variables.results_folder, 'calc_pwr.mat'));
            cfg_variables.MD_pwr=pwr_metrics.MD_pwr;
            cfg_variables.R_pwr=pwr_metrics.R_pwr;
            cfg_variables.op_per_joule=pwr_metrics.op_per_joule; 
        else
            cfg_variables.MD_pwr=missing;
            cfg_variables.R_pwr=missing;
            cfg_variables.op_per_joule=missing; 
        end

        if ~isfield(cfg_variables,'op_per_joule')
            cfg_variables.op_per_joule=missing;   
        end
    end
    
    cfg_variables_new=cfg_variables;
end