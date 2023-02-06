function [data_sim] = model_SET_TR(results_folder, netlist_folder, model_name, varargin)
    scale=1;
    kR='yes';
    imaxOpt=[1,1e-1,10];
    iminOpt=[3e-4,1e-5,5e-3];
    rsOpt=[25,1,200];
    etares_opt=50;
    vres_opt=-0.5;
    CH0_opt=1e-4;
    beta_opt=0.5;
    t0r_opt=1e4;
    v0r_opt=0.1;
    region='both';
    electroformed='no';
    r2sweep='RON';
    variability_model='no';
    scale_Roff=1e-3;

    
    
    prj_dir=fullfile('/home/users','aguirref','nn_rs_uab');
    root_folder=fullfile(prj_dir,'results','memdiode_sim');
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
            if strcmp(varargin{i},'etaresOpt')
                etaresOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'vresOpt')
                vresOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'etasetOpt')
                etasetOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'vsetOpt')
                vsetOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'betaOpt')
                betaOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'rsOpt')
                rsOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'rsmaxOpt')
                rsmaxOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'rsminOpt')
                rsminOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'alphamaxOpt')
                alphamaxOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'alphaminOpt')
                alphaminOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'region')
                region=varargin{i+1};
            end
            if strcmp(varargin{i},'r2sweep')
                r2sweep=varargin{i+1};
            end
            if strcmp(varargin{i},'rOFFscale')
                scale_Roff=varargin{i+1};
            end
            if strcmp(varargin{i},'electroformed')
                electroformed=varargin{i+1};
            end  
            if strcmp(varargin{i},'varModel')
                variability_model=varargin{i+1};
            end              
            if strcmp(varargin{i},'simTime')
                simulation_time=varargin{i+1};
            end           
            if strcmp(varargin{i},'voltages')
                voltages=varargin{i+1};
            end           
        end
    end   
    
    RUNLVL=5;

%     
%     figure();
%     ax1=subplot(1,2,1);
%     semilogy(Vt(:,2),abs(It(:,2)),'k');
%     ylabel('Memdiode current [A]');
%     xlabel('Memdiode Voltage [V]');
%     axis([-3 1.5 1e-7 1e-1]);
%     hold on
% 
%     ax2=subplot(1,2,2);
%     plot(Vt(:,2),abs(It(:,2)),'k');
%     ylabel('Memdiode current [A]');
%     xlabel('Memdiode Voltage [V]');
%     axis([-3 1.5 0 1.2e-2]);
%     hold on
multiWaitbar( 'CloseAll' );
multiWaitbar(sprintf('Model: %s...',model_name),0);

    for i=1:length(voltages)
        voltage=voltages(i);
        netlist_name='memdiode_SETtransient_netlist.sp';
        fileID = fopen(fullfile(netlist_folder,netlist_name),'w');
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

        fprintf(fileID,'.option nomod post newtol relmos=1e-7 absmos=1e-9\n');
        fprintf(fileID,'.model optmod opt itropt=300\n\n');

        fprintf(fileID,'.inc ''%s/models/memdiode/new_version/%s.sp''\n',prj_dir,model_name);
        fprintf(fileID,'Xmemdiodo pos gnd! memdiode H0=''H_init''\n');

        fprintf(fileID,'Rh h gnd! R=100k\n');
        fprintf(fileID,'.param H_init=0\n');

        fprintf(fileID,'Vsignal pos gnd! dc=0 pulse ( 0 %.3f 1u 100p 100p %.12e %.12e )\n',voltage, simulation_time,simulation_time);

        fprintf(fileID,'.PRINT TRAN V(pos) I(xmemdiodo.RS) V(h)\n');

        fprintf(fileID,'.MEAS TRAN final_val MAX I(xmemdiodo.RS)\n'); 
        fprintf(fileID,'.MEAS TRAN tr1 TRIG I(xmemdiodo.RS)=100u TD=0.2u RISE=1 TARG I(xmemdiodo.RS)=1m TD=0.2u RISE=1\n');
        fprintf(fileID,'.MEAS TRAN tset WHEN I(xmemdiodo.RS)=''final_val*0.1''\n');

        fprintf(fileID,'.tran 10u %f\n',simulation_time);

        fprintf(fileID,'.end\n');

        exe_file=strcat(hspice_root_dir,'bin',filesep,'hspice');
        system_command = strcat('module load HSPICE; "',exe_file,'" -i',{' '},fullfile(netlist_folder,netlist_name),{' '},'-n 1',{' '},'-o',{' '},fullfile(results_folder,strrep(netlist_name,'.sp','')));
        system_command=system_command{1};
        [status,cmdout]=unix(system_command,'-echo'); 

        
        sim_file=fopen(fullfile(results_folder,strrep(netlist_name,'.sp','.mt1.csv')));
        sim_file_line=fgetl(sim_file);
        while ~feof(sim_file)
            if ~isempty(strfind(sim_file_line,'final_val')) && ~isempty(strfind(sim_file_line,'tr1')) && ~isempty(strfind(sim_file_line,'tset'))
                headers=split(sim_file_line,',');
                sim_file_line=fgetl(sim_file);
                data=split(sim_file_line,',');
            else
                sim_file_line=fgetl(sim_file);
            end
        end
        
        for header_idx=1:length(headers)
            if strcmpi(headers{header_idx,1},'final_val')
                if isempty(strfind(data{header_idx,1},'failed'))
                    tran_data(i,1)=str2num(data{header_idx,1});
                else
                    tran_data(i,1)=0;
                end
            elseif strcmpi(headers{header_idx,1},'tr1')
                if isempty(strfind(data{header_idx,1},'failed'))
                    tran_data(i,2)=str2num(data{header_idx,1});
                else
                    tran_data(i,2)=0;
                end
            elseif strcmpi(headers{header_idx,1},'tset')
                if isempty(strfind(data{header_idx,1},'failed'))
                    tran_data(i,3)=str2num(data{header_idx,1});
                else
                    tran_data(i,3)=0;
                end
            end
        end
        fclose(sim_file);
        delete(fullfile(results_folder,strrep(netlist_name,'.sp','.mt1.csv')))
    
        data_sim{i,1}=importdata(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0')));
        
        multiWaitbar(sprintf('Model: %s...',model_name),i/length(voltages));
        
    end

    figure()
    subplot(2,2,1)
    semilogy(voltages,tran_data(:,1));
    xlabel('voltage [V]');
    ylabel('Final Current [A]');
    
    subplot(2,2,2)
    semilogy(voltages,tran_data(:,2));
    xlabel('voltage [V]');
    ylabel('Rise time [Sec]');
    
    subplot(2,2,3)
    semilogy(voltages,tran_data(:,3));
    xlabel('voltage [V]');
    ylabel('Time-to-set [A]');
    
    subplot(2,2,4)
    semilogy(voltages,1e-3./tran_data(:,3));
    xlabel('voltage [V]');
    ylabel('Transition Rate [A/sec]');

end