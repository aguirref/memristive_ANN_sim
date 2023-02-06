function [] = read_FSDB(results_folder, simEngine, wv_dir, netlist_name, neural_layers, sim_WR, remove_RAW_after_sim)
%IMPORTFSDB Summary of this function goes here
%   Detailed explanation goes here
    wv_commands=fopen(fullfile(results_folder,'read_fsdb.tcl'),'w');
    abs_path=what(results_folder);
    if strcmp(simEngine,'FineSim')
        fprintf(wv_commands,'set 	fileID 			[sx_open_sim_file_read "%s"]\n',fullfile(abs_path.path,strrep(netlist_name,'.net','.fsdb')));
    elseif strcmp(simEngine,'HSPICE')
        fprintf(wv_commands,'set 	fileID 			[sx_open_sim_file_read "%s"]\n',fullfile(abs_path.path,strrep(netlist_name,'.net','.tr1')));
    end
    
    if sim_WR==1
        fprintf(wv_commands,'set 	signals 		[ sx_signal *write_ok* *tim* ]\n');
    else
        for CPA_i=1:length(neural_layers)-1
            fprintf(wv_commands,'set 	signals 		[ sx_signal *neuron%s* *tim* ]\n',lower(neural_layers{CPA_i+1}));
            
%             fprintf(wv_commands,'foreach signal $signals {\n');
%             fprintf(wv_commands,'   set signal_name [ sx_signal_attribute $signal -name]\n');
%             fprintf(wv_commands,'   set file2save   "%s/$signal_name.m"\n',fullfile(abs_path.path,'signals'));
%             fprintf(wv_commands,'   sx_export_mfile $file2save $signal_name\n');
%             fprintf(wv_commands,'}\n');
%
%             Sadly this does not work: All signals are exported with a
%             different time vector...
        end
    end
    fprintf(wv_commands,'sx_export_mfile "%s" $signals',fullfile(abs_path.path,strrep(netlist_name,'.net','.m')));
    fclose(wv_commands);
    exe_file=strcat(wv_dir,'bin',filesep,'wv');
    mkdir(fullfile(abs_path.path,'signals'));

    system_command = strcat('',exe_file,' -ace_no_gui',{' '},fullfile(abs_path.path,'read_fsdb.tcl'));
    [status,cmdout]=system(system_command{1,1},'-echo');
    
    fileID=fopen(fullfile(abs_path.path,strrep(netlist_name,'.net','.m')),'r');
    while ~feof(fileID)
        current_line = fgetl(fileID);
        if ~contains(current_line,'Custom WaveView')
            pos_idx=strfind(current_line,' = ');
            signal_name=current_line(1:pos_idx-1);

            fileID_signal=fopen(fullfile(abs_path.path,'signals',sprintf('%s_file.m',signal_name)),'w');
            fprintf(fileID_signal,'%s',current_line);
            fclose(fileID_signal);
        end
    end
    fclose(fileID);

    delete(fullfile(results_folder,strrep(netlist_name,'.net','.m')));

    if strcmpi(remove_RAW_after_sim,'yes')
        if strcmpi(simEngine,'FineSim')
            delete(fullfile(results_folder,strrep(netlist_name,'.net','.fsdb')));   
        elseif strcmpi(simEngine,'HSPICE')
            delete(fullfile(results_folder,strrep(netlist_name,'.net','.tr1')));   
        end
    end
end

