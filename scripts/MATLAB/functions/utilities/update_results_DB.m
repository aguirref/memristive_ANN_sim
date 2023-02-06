function [] = update_results_DB(prov_results_directory,varargin)
    verify_settings=0;
    if mod(nargin,1)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-1
            if strcmp(varargin{i},'config')
                cfg_variables=varargin{i+1};
            end
            if strcmp(varargin{i},'results_directory')
                table_directory=varargin{i+1};
                if ispc
                    table_directory = strsplit(table_directory,'\');
                    table_directory = join(b(1:end-3),'\');  
                    table_directory = sprintf('%s\results\ANN_SPICE_sim');
                elseif isunix
                    table_directory = strsplit(table_directory,'/');
                    table_directory = join(b(1:end-3),'/');  
                    table_directory = sprintf('%s/results/ANN_SPICE_sim');
                end
                if ~exist(table_directory,'dir')
                    mkdir(table_directory);
                end
            end
            if strcmp(varargin{i},'verifySettings')
                verify_settings=varargin{i+1};
            end
        end
    end

    cfg_variables_original=cfg_variables;
    
    results_directory=fullfile(cfg_variables.project_directory,'results');

    filelist = dir(prov_results_directory);
    filelist = filelist(~[filelist.isdir]);

    j=1;
    for i=1:length(filelist)
        if strcmpi(filelist(i,1).name,'inference_results.mat')
            results_list(j,1)=filelist(i,1);
            j=j+1;
        end
    end

    if ~exist('results_list','var')
        filelist = dir(fullfile(results_directory, '**', 'MC=*','*.*'));  %get list of files and folders in any subfolder
        filelist = filelist(~[filelist.isdir]);  %remove folders from list
        
        j=1;
        for i=1:length(filelist)
            if strcmpi(filelist(i,1).name,'inference_results.mat')
                results_list(j,1)=filelist(i,1);
                j=j+1;
            end
        end
    end
   

    for i=1:length(results_list)

        load(fullfile(results_list(i,1).folder, results_list(i,1).name));
        load(fullfile(results_list(i,1).folder, 'inference_status.mat'));

        if exist(fullfile(results_list(i,1).folder, 'cfg_variables.mat'))
            load(fullfile(results_list(i,1).folder, 'cfg_variables.mat'));
 
            if verify_settings==1
                cfg_variables=update_cfg_variables(cfg_variables,results_list(i,1).folder);
                save(fullfile(results_list(i,1).folder,'cfg_variables.mat'),'cfg_variables');
            end
            if ~exist('table_directory','var')
                save2table(cfg_variables,results_directory);
            else
                save2table(cfg_variables,table_directory);
            end
        else
            if verify_settings==1
                cfg_variables=cfg_variables_original;
                cfg_variables.results_folder=results_list(i,1).folder;
                cfg_variables=update_cfg_variables(cfg_variables,results_list(i,1).folder);
            else
                cfg_variables.results_folder=results_list(i,1).folder;
                cfg_variables=update_cfg_variables(cfg_variables,results_list(i,1).folder,'mode','clean_only');
            end
            
            if strcmpi(cfg_variables.remap_W,'no')
                cfg_variables.remap_mode=[];
            end
            
            for hidden_layer_i=1:6
                if length(cfg_variables.neural_layers)<hidden_layer_i+2
                    eval(sprintf('cfg_variables.hidden_layer_%d=missing;',hidden_layer_i));
                end
            end

            if strcmpi(cfg_variables.train_tool,'Python')
                cfg_variables.mode=[];
                cfg_variables.numSigma=missing;
            end
            
            if ~iscell(cfg_variables.neural_layers{1,1})
                cfg_variables.neural_layers={cfg_variables.neural_layers};
            end
            if ~iscell(cfg_variables.partitions_N)
                cfg_variables.partitions_N={cfg_variables.partitions_N};
            elseif iscell(cfg_variables.partitions_N{1,1})
                cfg_variables.partitions_N=cfg_variables.partitions_N{1,1};
            end
            if ~iscell(cfg_variables.sys_clk)
                cfg_variables.sys_clk={cfg_variables.sys_clk};
            elseif iscell(cfg_variables.sys_clk{1,1})
                cfg_variables.sys_clk=cfg_variables.sys_clk{1,1};
            end
            if ~iscell(cfg_variables.sys_clk_WR_margin)
                cfg_variables.sys_clk_WR_margin={cfg_variables.sys_clk_WR_margin};
            elseif iscell(cfg_variables.sys_clk_WR_margin{1,1})
                cfg_variables.sys_clk_WR_margin=cfg_variables.sys_clk_WR_margin{1,1};
            end            
                       
            cfg_variables.accuracy_sim=1-SIM_DATA.ERROR_gral;
            cfg_variables.confusion_mat_sim=SIM_DATA.PROB;
            cfg_variables.avg_currents=SIM_DATA.CURR_MAT;
            cfg_variables.precision=SIM_DATA.PRECISION_gral;
            cfg_variables.sensitivity=SIM_DATA.SENSITIVITY_gral;
            cfg_variables.specificity=SIM_DATA.SPECIFICITY_gral;
            cfg_variables.f1_score=SIM_DATA.F1_score_gral;
    
            cfg_variables.sim_time=inference_vars.sim_time;
            cfg_variables.spice_sim_time=inference_vars.SPICE_sim_tim;
            cfg_variables.w_pos_init=inference_vars.W_pos_init;
            cfg_variables.w_neg_init=inference_vars.W_neg_init;
            cfg_variables.w_pos_init_fresh=inference_vars.W_pos_init_fresh;
            cfg_variables.w_neg_init_fresh=inference_vars.W_neg_init_fresh;
    
            if exist(fullfile(results_list(i,1).folder, 'calc_pwr.mat'),'file')
                pwr_metrics=load(fullfile(results_list(i,1).folder, 'calc_pwr.mat'));
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

            if strcmpi(cfg_variables.noise_transient,'yes')
                filelist_noise = dir(fullfile(cfg_variables.results_folder,'noise', '**', 'inference_results.mat'));  %get list of files and folders in any subfolder
                filelist_noise = filelist_noise(~[filelist_noise.isdir]);  %remove folders from list
                if ~isempty(filelist_noise)
                    for noise_folder_i=1:length(filelist_noise)
                        
                        noise_result_info=split(filelist_noise(noise_folder_i,1).folder,'/');
                        for j=1:length(noise_result_info)
                            str=noise_result_info{j,1};
                            if contains(str,'noise') && contains(str,'fmin') && contains(str,'fmax') && contains(str,'scale')
                                values=sscanf(str,'noise_fmax_%f_fmin_%f_scale_%f');
                                cfg_variables.noisefmax=values(1,1);
                                cfg_variables.noisefmin=values(2,1);
                                cfg_variables.noisescale=values(3,1);
                            end
                        end

                        load(fullfile(filelist_noise(noise_folder_i,1).folder, filelist_noise(noise_folder_i,1).name));
    
                        cfg_variables.SNR=SIM_DATA_noise.SNR;
                        cfg_variables.signal_rms=SIM_DATA_noise.SIGNAL_RMS;
                        cfg_variables.noise_rms=SIM_DATA_noise.NOISE_RMS;
                        cfg_variables.accuracy_noise=1-SIM_DATA_noise.ERROR_gral;
    
                        if ~exist('table_directory','var')
                            save2table(cfg_variables,results_directory);
                        else
                            save2table(cfg_variables,table_directory);
                        end
                    end
                else
                    cfg_variables.SNR=missing;
                    cfg_variables.signal_rms=missing;
                    cfg_variables.noise_rms=missing;
                    cfg_variables.accuracy_noise=missing;  
                    cfg_variables.noisefmax=missing;
                    cfg_variables.noisefmin=missing;
                    cfg_variables.noisescale=missing;
                                
                    if ~exist('table_directory','var')
                        save2table(cfg_variables,results_directory);
                    else
                        save2table(cfg_variables,table_directory);
                    end
                end
            else
                cfg_variables.SNR=missing;
                cfg_variables.signal_rms=missing;
                cfg_variables.noise_rms=missing;
                cfg_variables.accuracy_noise=missing;  
                cfg_variables.noisefmax=missing;
                cfg_variables.noisefmin=missing;
                cfg_variables.noisescale=missing;

                if ~exist('table_directory','var')
                    save2table(cfg_variables,results_directory);
                else
                    save2table(cfg_variables,table_directory);
                end
            end
            save(fullfile(results_list(i,1).folder,'cfg_variables.mat'),'cfg_variables');
        end
        clear cfg_variables
    end
end