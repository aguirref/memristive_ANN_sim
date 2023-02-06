function [return_Data] = sim_memdiode_dynamic(memdiode_model, model_ver, read_voltage, V_prog, tdelay_gral, tON_verify, wrFreq, tr, tOFF, tON_write, R_serie, target_curr, stimulus, varargin)
    plot_QMM=0;
    prj_dir=fullfile('/home/users','aguirref','nn_rs_uab');
    root_folder=fullfile(prj_dir,'results','memdiode_sim');
    results_folder='/home/users/aguirref/nn_rs_uab/results/memdiode_sim';
    simTime=5e-2;
    
    if mod(nargin-13,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-13
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
            if strcmp(varargin{i},'rootFolder')
                root_folder=varargin{i+1};
            end
            if strcmp(varargin{i},'plot_QMM')
                plot_QMM=varargin{i+1};
            end
        end
    end   
    if strcmpi(model_ver,'new')
        results_folder=fullfile(root_folder,'new_version',filesep,memdiode_model);
        models_folder=fullfile(prj_dir,'models','memdiode','new_version');
    elseif strcmpi(model_ver,'old')
        results_folder=fullfile(root_folder,'old_version',filesep,memdiode_model);
        models_folder=fullfile(prj_dir,'models','memdiode','old_version');
    elseif strcmpi(model_ver,'old_SBSF')
        results_folder=fullfile(root_folder,'old_SBSF_version',filesep,memdiode_model);
        models_folder=fullfile(prj_dir,'models','memdiode','old_SBSF_version');
    end
    
    if ~exist(results_folder,'dir')
        mkdir(results_folder)
    end
    
    netlist_name=strcat(memdiode_model,'.sp');
 
    hspice_root_dir='/usr/eda_tools/synopsys/HSPICE/O-2018.09/hspice/';
    
    RUNLVL=5;
    
    currents_above_target=[];
    while isempty(currents_above_target)
        
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

        if strcmpi(model_ver,'new')
            fprintf(fileID,'.inc ''%s.sp''\n',fullfile(models_folder,memdiode_model));
            fprintf(fileID,'Rserie Vpulsedsignal pos R=%.5g\n',R_serie);
            fprintf(fileID,'Xmemdiodo pos gnd! memdiode H0=0\n');
        elseif strcmpi(model_ver,'old')
            fprintf(fileID,'.inc ''%s.sp''\n',fullfile(models_folder,memdiode_model));
            fprintf(fileID,'Rserie Vpulsedsignal pos R=%.5g\n',R_serie);
            fprintf(fileID,'Xmemdiodo pos gnd! memdiode H0=0\n');        
        elseif strcmpi(model_ver,'old_SBSF')
            fprintf(fileID,'.inc ''%s.sp''\n',fullfile(models_folder,memdiode_model));
            fprintf(fileID,'Rserie Vpulsedsignal pos R=%.5g\n',R_serie);
            fprintf(fileID,'Xmemdiodo pos gnd! memdiode H0=1e-5\n');        
        end
        %fprintf(fileID,'Rh h gnd! R=1\n');
        %fprintf(fileID,'Vsetting-1 Vpulsedsignal Vpulsedsignal-1 dc=0 pulse ( 0 %f %.5g %.5g %.5g %.5g %.5g )\n',read_voltage, tdelay_gral, tr, tr, tON_verify, 1/wrFreq);    
        
        if strcmpi(stimulus,'pulses')
            fprintf(fileID,'Vsetting-2 Vpulsedsignal gnd! dc=0 pulse ( 0 %f %.5g %.5g %.5g %.5g %.5g )\n',V_prog, tdelay_gral+tOFF/2+tON_verify+tr, tr, tr, tON_write, 1/wrFreq);   
        elseif strcmpi(stimulus,'sin')
            fprintf(fileID,'Vsetting-2 Vpulsedsignal gnd! dc=0 sin ( 0 %.5g %.5g )\n', V_prog, wrFreq);   
        elseif strcmpi(stimulus,'ramp')
            pwl_signal=[0 0];
            
            while pwl_signal(end,1)<1/wrFreq
                if pwl_signal(end,1)<1/wrFreq/4
                    pwl_signal(end+1,:)=[1/wrFreq*1e-3+pwl_signal(end,1) 4*V_prog*1e-3+pwl_signal(end,2)]; 
                elseif pwl_signal(end,1)>1/wrFreq/4 && pwl_signal(end,1)<3/wrFreq/4
                    pwl_signal(end+1,:)=[1/wrFreq*1e-3+pwl_signal(end,1) -4*V_prog*1e-3+pwl_signal(end,2)]; 
                else
                    pwl_signal(end+1,:)=[1/wrFreq*1e-3+pwl_signal(end,1) 4*V_prog*1e-3+pwl_signal(end,2)]; 
                end
            end            
            fprintf(fileID,'Vsetting-2 Vpulsedsignal gnd! dc=0 PWL ( '); 
            for i=1:length(pwl_signal)
                fprintf(fileID,'%.4e %.4e ', pwl_signal(i,1), pwl_signal(i,2)); 
            end
            fprintf(fileID,')\n'); 
        end
        
        fprintf(fileID,'.PRINT TRAN V(pos) I(xmemdiodo.RS) V(Vpulsedsignal) V(Vpulsedsignal-1)\n');
        fprintf(fileID,'.MEAS TRAN measured_Imax param Xmemdiodo.Imax\n');
        fprintf(fileID,'.MEAS TRAN measured_Imin param Xmemdiodo.Imin\n');

        if strcmpi(stimulus,'pulses')
            fprintf(fileID,'.tran %.5g %.5g\n',tON_verify/5, simTime);
        elseif strcmpi(stimulus,'sin') || strcmpi(stimulus,'ramp')
            fprintf(fileID,'.tran %.5g %.5g\n',1/wrFreq*1e-3, 1/wrFreq);
        end
        

        fprintf(fileID,'.end\n');
        
        fclose(fileID);
        exe_file=strcat(hspice_root_dir,'bin',filesep,'hspice');
        system_command = strcat('module load HSPICE; "',exe_file,'" -i',{' '},fullfile(results_folder,netlist_name),{' '},'-n 1',{' '},'-o',{' '},fullfile(results_folder,strrep(netlist_name,'.sp','')));
        system_command=system_command{1};
        [status,cmdout]=unix(system_command,'-echo'); 

        repeat_data_extraction=1;
        while repeat_data_extraction==1
            data_sim=importdata(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0')));
            if exist('data_sim','var')
                if isfield(data_sim,'data')
                    if ~isempty(data_sim.data)
                         repeat_data_extraction=0;
                    end
                end
            end
        end
        
        sim_I=data_sim.data(:,3);
        sim_time=data_sim.data(:,1);
        sim_V=data_sim.data(:,2);
        %sim_Vpulsed=data_sim.data(:,4)-data_sim.data(:,5);

        if plot_QMM==1
            if strcmpi(stimulus,'sin') || strcmpi(stimulus,'ramp')
                figure('Name','memdiode_IV_characteristic');
                
                subplot(1,2,1)
                semilogy(sim_V,abs(sim_I),'-r');
                xlabel('Write time [sec]')
                ylabel('Memdiode Voltage [V]')
                hold on

                subplot(1,2,2)        
                plot(sim_V,sim_I,'-b');   
                %plot(sim_time,sim_I.*sim_Vpulsed,'-b');   
                xlabel('Write time [sec]')
                ylabel('Memdiode Current [A]')
                hold on               
            else

                subplot(1,2,1)
                plot(sim_time,sim_V,'-r');
                xlabel('Write time [sec]')
                ylabel('Memdiode Voltage [V]')
                hold on

                subplot(1,2,2)        
                plot(sim_time,sim_I,'-b');   
                %plot(sim_time,sim_I.*sim_Vpulsed,'-b');   
                xlabel('Write time [sec]')
                ylabel('Memdiode Current [A]')
                hold on
            end
        end

        if strcmpi(stimulus,'pulses')
            %[currents_above_target,~]=find((sim_I.*sim_Vpulsed)>target_curr);
            [currents_above_target,~]=find((sim_I)>target_curr);
            simTime=simTime*1.5;
            display('repeat!');
        else
            currents_above_target=1;
        end
    end
    
    if strcmpi(stimulus,'pulses')
        return_Data=sim_time(currents_above_target(1));
    else
        return_Data={sim_time, sim_V, sim_I};
    end
end