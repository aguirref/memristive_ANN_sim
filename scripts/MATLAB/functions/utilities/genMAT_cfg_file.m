function [cfg_vars] = genMAT_cfg_file(varargin)
    
    db_n_train_cfg=fullfile('..','config_files',filesep,'defaults',filesep,'DB_n_train_config.m');
    CPA_cfg=fullfile('..','config_files',filesep,'defaults',filesep,'CPA_config.m');
    conn_cfg=fullfile('..','config_files',filesep,'defaults',filesep,'conn_config.m');
    sgn_timming_cfg=fullfile('..','config_files',filesep,'defaults',filesep,'sgn_timing_config.m');
    sim_cfg=fullfile('..','config_files',filesep,'defaults',filesep,'sim_config.m');
    proc_settings=fullfile('..','config_files',filesep,'defaults',filesep,'process_settings.m');

    %% Optional values assignment
    if mod(nargin,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin
            if strcmp(varargin{i},'db_n_train_cfg')
                db_n_train_cfg=varargin{i+1};
            end
            if strcmp(varargin{i},'CPA_cfg')
               CPA_cfg=varargin{i+1};
            end
            if strcmp(varargin{i},'conn_cfg')
               conn_cfg=varargin{i+1};
            end    
            if strcmp(varargin{i},'sgn_timing_cfg')
               sgn_timming_cfg=varargin{i+1};
            end
            if strcmp(varargin{i},'sim_cfg')
               sim_cfg=varargin{i+1};
            end 
            if strcmp(varargin{i},'proc_settings')
               proc_settings=varargin{i+1};
            end 
        end
    end

    %% Load default program settings
    run(db_n_train_cfg);
    run(CPA_cfg);
    run(conn_cfg);
    run(sgn_timming_cfg);
    run(sim_cfg);    

    %% User-dependent system-configurable settings    
    run(proc_settings);  

    variables_list=who;

    for i=1:length(variables_list)
        if ~strcmpi(variables_list{i,1},'db_n_train_cfg') && ~strcmpi(variables_list{i,1},'CPA_cfg') && ~strcmpi(variables_list{i,1},'conn_cfg') && ~strcmpi(variables_list{i,1},'sgn_timming_cfg') && ~strcmpi(variables_list{i,1},'sim_cfg') && ~strcmpi(variables_list{i,1},'i')
            eval(sprintf('cfg_vars.%s=%s;',variables_list{i,1},variables_list{i,1}));
            eval(sprintf('clear %s',variables_list{i,1}));
        end
    end
end

