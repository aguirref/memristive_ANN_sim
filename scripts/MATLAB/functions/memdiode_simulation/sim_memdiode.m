function [I_limits] = sim_memdiode(memdiode_model, model_ver, read_voltage, varargin)
    scale=1;
    kR='yes';
    imaxOpt=[1,1e-1,10];
    iminOpt=[3e-4,1e-5,5e-3];
    rsOpt=[25,1,200];
    etares_opt='10';
    vres_opt='-0.5';
    CH0_opt='1e-3';
    beta_opt='0.5';
    t0r_opt='1e4';
    v0r_opt='0.1';
    region='SET';
    r2sweep='RON';
    plot_QMM=0;
    scale_Roff=1e-3;
    dev_polarity='positive';
    prj_dir=fullfile('/home/users','aguirref','nn_rs_uab');
    hspice_root_dir='/usr/eda_tools/synopsys/HSPICE/O-2018.09/hspice/';
    
    if mod(nargin-3,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-3
            if strcmp(varargin{i},'scale')
                scale=varargin{i+1};
            end
            if strcmp(varargin{i},'keepRatio')
                kR=varargin{i+1};
            end
            if strcmp(varargin{i},'imaxOpt')
                imaxOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'iminOpt')
                iminOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'rsOpt')
                rsOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'region')
                region=varargin{i+1};
            end
            if strcmp(varargin{i},'r2sweep')
                r2sweep=varargin{i+1};
            end
            if strcmp(varargin{i},'plot_QMM')
                plot_QMM=varargin{i+1};
            end
            if strcmp(varargin{i},'prj_dir')
                prj_dir=varargin{i+1};
            end
            if strcmp(varargin{i},'simEngine')
                simEngine=varargin{i+1};
            end
            if strcmp(varargin{i},'hspice_dir')
                hspice_root_dir=varargin{i+1};
            end
            if strcmp(varargin{i},'timeSTEP_min')
                timestep_sim=varargin{i+1};
            end
             if strcmp(varargin{i},'dev_polarity')
                dev_polarity=varargin{i+1};
            end
        end
    end   

    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
        model_syntax_ver='HSPICE';
    elseif strcmp(simEngine,'LTspice')
        model_syntax_ver='LTspice';
    else
        model_syntax_ver='HSPICE';
    end

    root_folder=fullfile(prj_dir,'results','memdiode_sim');

    if strcmpi(model_ver,'DMM')
        results_folder=fullfile(root_folder,'DMM',filesep,memdiode_model);
        models_folder=fullfile(prj_dir,'models','memdiode','DMM',model_syntax_ver);
    elseif strcmpi(model_ver,'QMM')
        results_folder=fullfile(root_folder,'QMM',filesep,memdiode_model);
        models_folder=fullfile(prj_dir,'models','memdiode','QMM',model_syntax_ver);
    elseif strcmpi(model_ver,'QMM_SBSF')
        results_folder=fullfile(root_folder,'QMM_SBSF',filesep,memdiode_model);
        models_folder=fullfile(prj_dir,'models','memdiode','QMM_SBSF',model_syntax_ver);
    end
    
    if ~exist(results_folder,'dir')
        mkdir(results_folder)
    end
    
    netlist_name=strcat(memdiode_model,'.sp');

    if exist(fullfile(prj_dir,'models','memdiode',model_ver,model_syntax_ver,memdiode_model),'dir')
        file_Vt=fullfile(prj_dir,'models','memdiode',model_ver,model_syntax_ver,memdiode_model,'mean_Vt_data.csv');
        file_Vt_reset=fullfile(prj_dir,'models','memdiode',model_ver,model_syntax_ver,memdiode_model,'mean_Vt_data.csv');    
        file_Vt_complete=fullfile(prj_dir,'models','memdiode',model_ver,model_syntax_ver,memdiode_model,'mean_Vt_data.csv');        
    else
        file_Vt=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_data_set.txt');
        file_Vt_reset=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_data_reset.txt');    
        file_Vt_complete=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_complete_2.txt');
    end

    if strcmpi(region,'SET') || strcmpi(region,'both')
        Vt=importdata(file_Vt);
    end
    
    if strcmpi(region,'RESET')
        Vt_reset=importdata(file_Vt_reset);
    end
    
    RUNLVL=5;
        
    if strcmpi(region,'SET') || strcmpi(region,'both')    
        fileID = fopen(fullfile(results_folder,netlist_name),'w');
        fprintf(fileID,'*******************************************************\n');

        fprintf(fileID,'.option INGOLD=1\n');
        fprintf(fileID,'.option RUNLVL=%d\n',RUNLVL);
        fprintf(fileID,'.option SPLIT_DP=1\n');
        fprintf(fileID,'.option MEASFORM=3\n');
        fprintf(fileID,'.option CSDF=1\n');    
        fprintf(fileID,'.option METHOD=TRAP\n');
        fprintf(fileID,'.option LIST=0\n');
        fprintf(fileID,'.option LIS_NEW=1\n'); 
        fprintf(fileID,'.option POST=CSDF\n');    

        connections = 'pos gnd!';
        if strcmpi(dev_polarity,'positive')
            multiplier_rVoltage=1;
        else
            multiplier_rVoltage=-1;
        end

        if strcmpi(model_ver,'DMM')
            fprintf(fileID,'.inc ''%s.sp''\n',fullfile(models_folder,memdiode_model));
            fprintf(fileID,'Xmemdiodo %s memdiode H0=0\n', connections);
        elseif strcmpi(model_ver,'QMM')
            fprintf(fileID,'.inc ''%s.sp''\n',fullfile(models_folder,memdiode_model));
            fprintf(fileID,'Xmemdiodo %s memdiode H0=0\n', connections);        
        elseif strcmpi(model_ver,'QMM_SBSF')
            fprintf(fileID,'.inc ''%s.sp''\n',fullfile(models_folder,memdiode_model));
            fprintf(fileID,'Xmemdiodo %s memdiode H0=1e-5\n', connections);        
        end
        %fprintf(fileID,'Rh h gnd! R=1\n');
        fprintf(fileID,'Vsignal pos gnd! PWL PWLFILE=''%s''\n',file_Vt_complete);

        fprintf(fileID,'.PRINT TRAN V(pos) I(xmemdiodo.RS) V(xmemdiodo.H)\n');
        fprintf(fileID,'.MEAS TRAN measured_Imax param Xmemdiodo.Imax\n');
        fprintf(fileID,'.MEAS TRAN measured_Imin param Xmemdiodo.Imin\n');

        fprintf(fileID,'.tran %.3g %f\n',timestep_sim,Vt(end,1));

        fprintf(fileID,'.end\n');

        exe_file=strcat(hspice_root_dir,'bin',filesep,'hspice');
        system_command = strcat('module load HSPICE; "',exe_file,'" -i',{' '},fullfile(results_folder,netlist_name),{' '},'-n 1',{' '},'-o',{' '},fullfile(results_folder,strrep(netlist_name,'.sp','')));
        system_command=system_command{1};
        [status,cmdout]=unix(system_command,'-echo'); 

        if exist(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr1')),'file')
            data_sim=importdata(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr1')));
        elseif exist(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0')),'file')
            data_sim=importdata(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0')));
        end

        sim_I=data_sim.data(:,3);
        sim_time=data_sim.data(:,1);
        sim_V=data_sim.data(:,2);
        sim_H=data_sim.data(:,4);

        if plot_QMM==1
        
            figure('Name','memdiode_IV_characteristic');

            subplot(1,2,1)
            semilogy(sim_V,abs(sim_I),'--r');
            xlabel('Memdiode Voltage [V]')
            ylabel('Memdiode Current [A]')
            hold on

            subplot(1,2,2)        
            plot(sim_V,sim_H,'--r');   
            xlabel('Memdiode Voltage [V]')
            ylabel('Memdiode Memory [a.u.]')
            hold on
        end
        
        [~,aux_var_pos]=findpeaks(sim_V*multiplier_rVoltage);

        sim_V_HRS=sim_V(1:aux_var_pos);
        sim_I_HRS=sim_I(1:aux_var_pos);
        aux_pos_minor=find(sim_V_HRS<read_voltage*multiplier_rVoltage);
        aux_pos_major=find(sim_V_HRS>read_voltage*multiplier_rVoltage);
        if strcmpi(dev_polarity,'positive')
            IMIN=interp1(sim_V_HRS([max(aux_pos_minor) min(aux_pos_major)]),sim_I_HRS([max(aux_pos_minor) min(aux_pos_major)]),read_voltage*multiplier_rVoltage);
        else
            IMAX=abs(interp1(sim_V_HRS([max(aux_pos_major) min(aux_pos_minor)]),sim_I_HRS([max(aux_pos_major) min(aux_pos_minor)]),read_voltage*multiplier_rVoltage));
        end

        sim_V_LRS=sim_V(aux_var_pos:end);
        sim_I_LRS=sim_I(aux_var_pos:end);
        aux_pos_minor=find(sim_V_LRS>read_voltage*multiplier_rVoltage);
        aux_pos_major=find(sim_V_LRS<read_voltage*multiplier_rVoltage);
        if strcmpi(dev_polarity,'positive')
            IMAX=interp1(sim_V_LRS([max(aux_pos_minor) min(aux_pos_major)]),sim_I_LRS([max(aux_pos_minor) min(aux_pos_major)]),read_voltage*multiplier_rVoltage);
        else
            IMIN=abs(interp1(sim_V_LRS([max(aux_pos_major) min(aux_pos_minor)]),sim_I_LRS([max(aux_pos_major) min(aux_pos_minor)]),read_voltage*multiplier_rVoltage));
        end

        I_limits=[IMIN*1.02 IMAX*0.98];
    end
    

end