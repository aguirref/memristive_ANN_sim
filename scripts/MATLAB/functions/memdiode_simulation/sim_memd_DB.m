function [sim_data] = sim_memd_DB(beta,...
                                  etaset,...
                                  vset,...
                                  etares,...
                                  vres,...
                                  imax,...
                                  alpha,...
                                  rs,...
                                  imin,...
                                  RUNLVL,...
                                  tDELMAX,...
                                  Vapp,...
                                  tstep,...
                                  tstop,...
                                  re_freq,...
                                  vt_mode,...
                                  varargin)

    results_folder='/home/users/aguirref/nn_rs_uab/results/memdiode_sim';
    netlist_folder='/home/users/aguirref/nn_rs_uab/netlists';
    hspice_root_dir='/usr/eda_tools/synopsys/HSPICE/O-2018.09/hspice/';
    prj_dir=fullfile('/home/users','aguirref','nn_rs_uab');
    root_folder=fullfile(prj_dir,'results','memdiode_sim');
    useAlphaMaxMin=0;
    file_Vt_complete=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_complete.txt');
    
    if mod(nargin-16,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-16
            if strcmp(varargin{i},'alphamax')
                alphamax=varargin{i+1};
                useAlphaMaxMin=1;
            end
            if strcmp(varargin{i},'alphamin')
                alphamin=varargin{i+1};
                useAlphaMaxMin=1;
            end        
        end
    end   
    
   
    netlist_name='memdiode_fit_netlist.sp';
    fileID = fopen(fullfile(netlist_folder,netlist_name),'w');
    fprintf(fileID,'*******************************************************\n');

    fprintf(fileID,'.option INGOLD=1\n');
    fprintf(fileID,'.option RUNLVL=%d\n',RUNLVL);
    fprintf(fileID,'.option DELMAX=%.2e\n',tDELMAX);
    fprintf(fileID,'.option SPLIT_DP=1\n');
    fprintf(fileID,'.option MEASFORM=3\n');
    fprintf(fileID,'.option CSDF=1\n');    
    fprintf(fileID,'.option METHOD=TRAP\n');
    fprintf(fileID,'.option LIST=0\n');
    fprintf(fileID,'.option LIS_NEW=1\n'); 
    fprintf(fileID,'.option POST=CSDF\n');    

    if strcmpi(version,'new')
        fprintf(fileID,'.inc ''%s/models/memdiode/model_fit/memdiode_HSPICE_new_base.sp''\n',prj_dir);
        fprintf(fileID,'Xmemdiodo pos gnd! h memdiode H0=''H_init'' beta=''beta_opt'' T0s=''T0s_opt'' V0s=''V0s_opt'' T0r=''T0r_opt'' V0r=''V0r_opt'' imax=''imax_opt'' alphamax=''alphamax_opt'' rsmax=''rs_opt'' imin=''imin_opt'' alphamin=''alphamin_opt'' rsmin=''rs_opt''\n');
    elseif strcmpi(version,'old_SBSF')
        fprintf(fileID,'.inc ''%s/models/memdiode/model_fit/memdiode_HSPICE_old_SBSF_base.sp''\n',prj_dir);
        fprintf(fileID,'Xmemdiodo pos gnd! h memdiode H0=''H_init'' H0_limit=%.3e beta=''beta_opt'' etaset=''etaset_opt'' vset=''vset_opt'' vt=''vt_opt'' etares=''etares_opt'' vres=''vres_opt'' imax=''imax_opt'' alphamax=''alphamax_opt'' rsmax=''rsmax_opt'' imin=''imin_opt'' alphamin=''alphamin_opt'' rsmin=''rsmin_opt'' isb=''isb_opt'' gam=''gam_opt''\n', H0_limit);        
    else
        fprintf(fileID,'.inc ''%s/models/memdiode/model_fit/memdiode_HSPICE_old_base.sp''\n',prj_dir);
        if useAlphaMaxMin==0
            fprintf(fileID,'Xmemdiodo pos gnd! h memdiode H0=0 beta=%.1e etaset=%.1e vset=%.1e etares=%.1e vres=%.1e imax=%.1e alphamax=%.1e rsmax=%.1e imin=%.1e alphamin=%.1e rsmin=%.1e CH0=%.1e\n',beta, etaset, vset, etares, vres, imax, alpha, rs, imin, alpha, rs, 1e-3);        
        else
            fprintf(fileID,'Xmemdiodo pos gnd! h memdiode H0=0 beta=%.1e etaset=%.1e vset=%.1e etares=%.1e vres=%.1e imax=%.1e alphamax=%.1e rsmax=%.1e imin=%.1e alphamin=%.1e rsmin=%.1e CH0=%.1e\n',beta, etaset, vset, etares, vres, imax, alphamax, rs, imin, alphamin, rs, 1e-3);                    
        end
    end
    fprintf(fileID,'Rh h gnd! R=100k\n');

    if strcmpi(vt_mode,'sin')
        fprintf(fileID,'Vsignal pos gnd! SIN(0 %.2e)\n',Vapp);
    elseif strcmpi(vt_mode,'vt_data')
        fprintf(fileID,'Vsignal pos gnd! PWL PWLFILE=''%s''\n',file_Vt_complete);
    end
    
    fprintf(fileID,'.PRINT TRAN V(pos) I(xmemdiodo.RS) V(h)\n');
    fprintf(fileID,'.tran %.2e %.2e\n',tstep, tstop);
    fprintf(fileID,'.end\n');
    fclose(fileID);
    
    exe_file=strcat(hspice_root_dir,'bin',filesep,'hspice');
    system_command = strcat('module load HSPICE; "',exe_file,'" -i',{' '},fullfile(netlist_folder,netlist_name),{' '},'-n 1',{' '},'-o',{' '},fullfile(results_folder,strrep(netlist_name,'.sp','')));
    system_command=system_command{1};
    [status,cmdout]=unix(system_command,'-echo'); 

    fid = fopen(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0')));
    fid2 = fopen(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0-pre')),'w+');
    fid3 = fopen(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0-post')),'w+');

    tline = fgetl(fid);
    lineCounter = 1;
    first_line=999e3;
    last_line=999e3;
    printtri=1;
    while ischar(tline)
        %disp(tline)
        if strcmp(tline, 'x')
            first_line=lineCounter;
            last_line=999e3;
            %break;
        end

        if lineCounter>first_line && lineCounter<last_line
            if printtri==1
                fprintf(fid2,'%s\n',tline);
            else
                fprintf(fid3,'%s\n',tline);
            end
        end

        tline = fgetl(fid);
        lineCounter = lineCounter + 1; 

        if strcmp(tline, 'y')
            last_line=lineCounter;
            first_line=999e3;
            printtri=printtri+1;
            %break;
        end       
    end
    fclose(fid);
    fclose(fid2);
    fclose(fid3);

    data_sim_pre=importdata(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0-pre')));

    IV_sim_data=data_sim_pre.data;

    [v_app,time] = resample(IV_sim_data(:,2),IV_sim_data(:,1),re_freq,1,1);
    [I_memd,~] = resample(IV_sim_data(:,3),IV_sim_data(:,1),re_freq,1,1);
    [v_h,~] = resample(IV_sim_data(:,4),IV_sim_data(:,1),re_freq,1,1);
    
%     figure();
%     semilogy(v_app,abs(I_memd));
    
    sim_data=horzcat(time,v_app,I_memd,v_h);
    end