function map=simRRAM_MLP_dualPart(neu_per_layer,varargin)
    
    %% Default values    
    % Neural Network characteristics
    neural_layers={'A' 'B'};
    pos_neg={'pos' 'neg'};
    simulation_time=150e-3;
    max_time_step=0;

    %CPA connectivity
    connections='dualSide_connect_inputs';
    partitions_N=4;
    V_prog=1.25;
    VDD=1.2;
    read_voltage=0.5;
    wrFreq=100e3;
    tr=100e-9;
    DC=0.5;
    DC_v=0.5;
    tdelay_gral=10e-6;
    sel_type='ideal'; 
    dev_polarity = 'positive';

    % Parasitic resistances
    R_shunt=10;
    Rs=[1]; 
    R_cs=1;
    use_caps='no';
    C_memdiode=45e-18;
    C_interline=21.6e-15;
    Gsense='yes';
    TIA='no';
    R_sense=10e3;
    sigmaImax=0;
    sigmaImin=0;
    print_V_tolerance=1e-4;
    print_I_tolerance=1e-12;

    %simulation options    
    calc_pwr='no';
    calc_WR_margin=0;
    calc_RE_latency=0;
    genPlot='yes';
    RUNLVL=3;
    simulate_netlist='yes';
    runMode='-b';
    sim_WR=0;
    memdiode_model='memdiode3';
    simEngine='FineSim';
    noise_transient='no';
    max_pulses=150;
    
    hspice_root_dir='/usr/eda_tools/synopsys/HSPICE/O-2018.09/hspice/';
    if isunix
        sim_folder='/home/Simulations/local/Neural_network';
    else
        sim_folder='D:\NN_RS_UAB\Simulations';
    end
    prj_dir=fullfile('/home/users','aguirref','nn_rs_uab');    
    nProc=1;
    use_result_folder=0;
    input_vector_freq=100e3;
    in_vect_r_f_time=100e-9;
    inputMatrix=ones(neu_per_layer(1),1)*read_voltage;
    meas_statement='.MEAS';
    analysis='transient';
    read_delay=1e-6;
    correct_folder='no_correction';

    %CPA initialization
    W_matrix_pos=[];
    W_matrix_neg=[];
    W_matrix_init_pos=[];
    W_matrix_init_neg=[];
    guiStatus_data=0;
    current_dir=pwd;
    error_band=0.01;
    
    % save defaults
    save_write_signals='no';
    save_h='no';
    save_vpn='no';
    save_inner_neuron='no';
    save_imemd='no';
    debug_mode='no';
    save_all=0;
    output_fmt='ascii';
    
    if ispc
        current_dir=split(current_dir,"\");
        current_dir=current_dir(1:end-3);
        root = cell2mat(join(current_dir,"\"));
        root_str=strrep(root,'\','\\');
    elseif isunix
        current_dir=split(current_dir,"/");
        current_dir=current_dir(1:end-3);
        root_str = cell2mat(join(current_dir,"/"));
        %root_str=strrep(root,'\','\\');       
    end 
    
    %% Optional values assignment
    if mod(nargin-1,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-1
            if strcmp(varargin{i},'simTime')
                simulation_time=varargin{i+1};
                if simulation_time<10e-6
                    fprintf('---> Simulation time is too small and the circuit hasn''t event gone thorugh the initial reset. Chaging to 10e-6\n')
                    %simulation_time=10e-6;
                end
            end
            if strcmp(varargin{i},'hspice_dir')
                hspice_root_dir=varargin{i+1};
            end
            if strcmp(varargin{i},'finesim_dir')
                finesim_dir=varargin{i+1};
            end
            if strcmp(varargin{i},'wv_dir')
                wv_dir=varargin{i+1};
            end      
            if strcmp(varargin{i},'license_address')
                license_address=varargin{i+1};
                setenv('SNPSLMD_LICENSE_FILE',license_address);
                setenv('LM_LICENSE_FILE',license_address);                
            end               
            if strcmp(varargin{i},'minStep')
                max_time_step=varargin{i+1};
            end
            if strcmp(varargin{i},'Rs')
                Rs=varargin{i+1};
            end
            if strcmp(varargin{i},'model_ver')
                model_ver=varargin{i+1};
            end
            if strcmp(varargin{i},'Rcs')
                R_cs=varargin{i+1};
            end 
            if strcmp(varargin{i},'deg_folder')
                deg_folder=varargin{i+1};
            end
            if strcmp(varargin{i},'selType')
                sel_type=varargin{i+1};
            end            
            if strcmp(varargin{i},'CPAsettings')
                CPA_settings=varargin{i+1};
                neurons=CPA_settings.neurons;
                dev_polarity=CPA_settings.dev_polarity; 
            end            
            if strcmp(varargin{i},'SIMsettings')
                SIM_settings=varargin{i+1};
                save_write_signals=SIM_settings.save_write_signals;
                save_h=SIM_settings.save_h;
                save_vpn=SIM_settings.save_vpn;
                save_inner_neuron=SIM_settings.save_inner_neuron;
                save_imemd=SIM_settings.save_imemd;
                save_all=SIM_settings.save_all;
                debug_mode=SIM_settings.debug_mode;
                input_vector_freq=SIM_settings.input_vector_freq;
                in_vect_r_f_time=SIM_settings.in_vect_r_f_time;     
                noisefmax=SIM_settings.noisefmax;
                noisefmin=SIM_settings.noisefmin;
                noisescale=SIM_settings.noisescale;
                in_polarity=SIM_settings.in_polarity;
                wrFreq=SIM_settings.wrFreq;
                tr=SIM_settings.tr;
                DC=SIM_settings.DC;
                DC_v=SIM_settings.DC_v;
                tdelay_gral=SIM_settings.tdelay_gral;                
            end       
            if strcmp(varargin{i},'in_vect_r_f_time')
               in_vect_r_f_time=varargin{i+1};
            end
            if strcmp(varargin{i},'input_vector_freq')
               input_vector_freq=varargin{i+1};
            end
            if strcmp(varargin{i},'dir_n_files')
                dir_n_files=varargin{i+1};
            end
            if strcmp(varargin{i},'Rs_PAD')
                Rs_PAD=varargin{i+1};
            end 
            if strcmp(varargin{i},'R_shunt')
                R_shunt=varargin{i+1};
            end
            if strcmp(varargin{i},'calcWRmargin')
                calc_WR_margin=varargin{i+1};
            end            
            if strcmp(varargin{i},'calcRElatency')
                calc_RE_latency=varargin{i+1};
            end            
            if strcmp(varargin{i},'rVoltage') 
                read_voltage=varargin{i+1};
            end
            if strcmp(varargin{i},'sigmaImax')
                sigmaImax=varargin{i+1};
            end
            if strcmp(varargin{i},'sigmaImin')
                sigmaImin=varargin{i+1};
            end
            if strcmp(varargin{i},'W_matrix')
                W_matrix=varargin{i+1};
            end
            if strcmp(varargin{i},'calcPWR')
                calc_pwr=varargin{i+1};
            end
            if strcmp(varargin{i},'errorBand')
                error_band=varargin{i+1};
            end
            if strcmp(varargin{i},'removeRAW')
                remove_RAW_after_sim=varargin{i+1};
            end
            if strcmp(varargin{i},'nPart')
                partitions_N=varargin{i+1};
            end
            if strcmp(varargin{i},'W_matrix_pos')
                W_matrix_pos=varargin{i+1};
            end
            if strcmp(varargin{i},'W_matrix_neg')
                W_matrix_neg=varargin{i+1};
            end            
            if strcmp(varargin{i},'W_init')
                W_matrix_init=varargin{i+1};
            end 
            if strcmp(varargin{i},'W_init_pos')
                W_matrix_init_pos=varargin{i+1};
            end 
            if strcmp(varargin{i},'useWData')
                h_data=varargin{i+1};
            end
            if strcmp(varargin{i},'useCaps')
               use_caps=varargin{i+1};
            end
            if strcmp(varargin{i},'noiseTransient')
               noise_transient=varargin{i+1};
            end
            if strcmp(varargin{i},'cInterline')
               C_interline=varargin{i+1};
            end
            if strcmp(varargin{i},'senseMode')
                if strcmpi(varargin{i+1},'Gsense')
                    Gsense='yes';
                    TIA='no';
                elseif strcmpi(varargin{i+1},'TIA')
                    TIA='yes';
                    Gsense='no';
                end
            end
            if strcmp(varargin{i},'rSense')
               R_sense=varargin{i+1};
            end
            if strcmp(varargin{i},'tiaGain')
               TIA_gain=varargin{i+1};
            end
            if strcmp(varargin{i},'rinTIA')
               R_in_TIA=varargin{i+1};
            end
            if strcmp(varargin{i},'cMemdiode')
               C_memdiode=varargin{i+1};
            end
            if strcmp(varargin{i},'cLine2gnd')
               C_line2gnd=varargin{i+1};
            end
            if strcmp(varargin{i},'W_init_neg')
                W_matrix_init_neg=varargin{i+1};
            end 
            if strcmp(varargin{i},'simNetlist')
                simulate_netlist=varargin{i+1};
            end
            if strcmp(varargin{i},'runMode')
                runMode=varargin{i+1};
            end
            if strcmp(varargin{i},'layers')
                neural_layers=varargin{i+1};
            end            
            if strcmp(varargin{i},'guiStatus')
                guiStatus_data=~isempty(varargin{i+1});
                status_gui=varargin{i+1}.status_gui;
                str_to_display=varargin{i+1}.str_to_display;
                status_guiObject=varargin{i+1}.status_guiObject;
                str_header=varargin{i+1}.str_header;
            end
            if strcmp(varargin{i},'strDisplay')
                str_to_display=varargin{i+1};
            end
            if strcmp(varargin{i},'strHdr')
                str_header=varargin{i+1};
            end
            if strcmp(varargin{i},'simWR')
                sim_WR=varargin{i+1};
            end
            if strcmp(varargin{i},'memdiodeModel')
                memdiode_model=varargin{i+1};
            end
            if strcmp(varargin{i},'simEngine')
                simEngine=varargin{i+1};
            end  
            if strcmp(varargin{i},'analysis')
                analysis=varargin{i+1};
            end  
            if strcmp(varargin{i},'RUNLVL')
                RUNLVL=varargin{i+1};
            end              
            if strcmp(varargin{i},'genPlot')
                genPlot=varargin{i+1};
            end  
            if strcmp(varargin{i},'vProg')
                V_prog=varargin{i+1};
            end     
            if strcmp(varargin{i},'finesim_vprbtol')
                print_V_tolerance=varargin{i+1};
            end   
            if strcmp(varargin{i},'finesim_iprbtol')
                print_I_tolerance=varargin{i+1};
            end   
            if strcmp(varargin{i},'nProc')
                nProc=varargin{i+1};
            end  
            if strcmp(varargin{i},'sysClk')
                sys_clk=varargin{i+1};
            end  
            if strcmp(varargin{i},'saveALL')
                save_all=varargin{i+1};
            end 
            if strcmp(varargin{i},'connections')
                connections=varargin{i+1};
            end 
            if strcmp(varargin{i},'correctsRs')
                if strcmpi(varargin{i+1},'yes')
                    correct_folder='rs_corrected';
                else
                    correct_folder='no_correction';
                end
            end 
            if strcmp(varargin{i},'outputFMT')
                output_fmt=varargin{i+1};
            end  
            if strcmp(varargin{i},'inputMatrix')
                if size(varargin{i+1},1)==neu_per_layer(1)
                    inputMatrix=varargin{i+1};
                else
                    fprintf('The length of the input vector is different to the expected number of inputs. Using the default value at every input')
                    inputMatrix=ones(neu_per_layer(1),1)*read_voltage;
                end
            end
            if strcmp(varargin{i},'simDir')
                if varargin{i+1}==1
                    use_result_folder=1;
                elseif ischar(varargin{i+1})
                    sim_folder=varargin{i+1};
                end
            end  
            if strcmp(varargin{i},'results_folder')
                results_folder=varargin{i+1};
            end  
            Rs_PAD=Rs;
            neu_per_layer_sub=neu_per_layer./partitions_N;
        end
    end
    if use_result_folder==1
        sim_folder=fullfile(root_str,'results','Neural_Network',strcat('NW_',num2str(neu_per_layer(1)),'by',num2str(neu_per_layer(2)),'_closed_loop'),simEngine);
    elseif ~exist(sim_folder,'dir')
        mkdir(sim_folder);
    end
    if calc_RE_latency==1
        R_off_value='100G';
        unselected_row_voltage=V_prog/3;
        unselected_col_voltage=V_prog/3;        
    else
        R_off_value='100Meg';
        unselected_row_voltage=0;
        unselected_col_voltage=V_prog/2;
    end   
    %% Weights matrix generation
    if sim_WR==0
        if isempty(W_matrix_init_pos)
            W_matrix_init_pos=W_matrix_init;
            W_matrix_init_pos(W_matrix_init_pos<0)=0;
        end
        if isempty(W_matrix_init_neg)
            W_matrix_init_neg=W_matrix_init;
            W_matrix_init_neg(W_matrix_init_neg>0)=0;
            W_matrix_init_neg=W_matrix_init_neg*(-1);
            W_matrix_init(:,:,2)=W_matrix_init_neg;
        end
        for aux_index_ii=1:length(W_matrix_init_pos)
            W_matrix_init_aux(:,:,1)=W_matrix_init_pos{aux_index_ii,1};
            W_matrix_init_aux(:,:,2)=W_matrix_init_neg{aux_index_ii,1};
            W_matrix_init{1,aux_index_ii}=W_matrix_init_aux;    
            clear W_matrix_init_aux
            
            W_matrix_aux(:,:,1)=zeros(size(W_matrix_init_pos{aux_index_ii,1}));
            W_matrix_aux(:,:,2)=zeros(size(W_matrix_init_pos{aux_index_ii,1}));
            W_matrix{1,aux_index_ii}=W_matrix_aux; 
            clear W_matrix_aux            
        end
    else
        if isempty(W_matrix_pos)
            W_matrix_pos=W_matrix;
            W_matrix_pos(W_matrix_pos<=0)=10e-6;
        end
        if isempty(W_matrix_neg)
            W_matrix_neg=W_matrix;
            W_matrix_neg(W_matrix_neg>=0)=-10e-6;
            W_matrix_neg=W_matrix_neg*(-1);   
        end    
        for aux_index_ii=1:length(W_matrix_pos)
            W_matrix_aux(:,:,1)=W_matrix_pos{aux_index_ii,1};
            W_matrix_aux(:,:,2)=W_matrix_neg{aux_index_ii,1};
            W_matrix{1,aux_index_ii}=W_matrix_aux;    
            clear W_matrix_aux
            
            W_matrix_init_aux(:,:,1)=zeros(size(W_matrix_pos{aux_index_ii,1}));
            W_matrix_init_aux(:,:,2)=zeros(size(W_matrix_pos{aux_index_ii,1}));
            W_matrix_init{1,aux_index_ii}=W_matrix_init_aux; 
            clear W_matrix_init_aux         
            
        end
    end
    
    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
        if sim_WR==1 && neu_per_layer(1)*neu_per_layer(2)>225 && save_all==1 && strcmp(output_fmt,'ascii')
            fprintf('---> WARNING! The number of elements in the array is too large for printing the transient results to a .txt file. Using a binary format\n')
            output_fmt='binary';
        end
    end
    if simulation_time<1/input_vector_freq*size(inputMatrix,2)
        fprintf('---> WARNING! The simulation time is smaller than the time required to test all the input vectors considered. Consder using a simulation time of at least %.12e',1/input_vector_freq*size(inputMatrix,2));
    end
    
    %% Model parameters extraction
        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
            model_syntax_ver='HSPICE';
        elseif strcmp(simEngine,'LTspice')
            model_syntax_ver='LTspice';
        else
            model_syntax_ver='HSPICE';
        end

        if strcmpi(model_ver,'DMM')
            models_folder=fullfile(prj_dir,'models','memdiode','DMM',model_syntax_ver);
        elseif strcmpi(model_ver,'QMM')
            models_folder=fullfile(prj_dir,'models','memdiode','QMM',model_syntax_ver);
        elseif strcmpi(model_ver,'QMM_SBSF')
            models_folder=fullfile(prj_dir,'models','memdiode','QMM_SBSF',model_syntax_ver);    
        elseif strcmpi(model_ver,'yakopic')
            models_folder=fullfile(prj_dir,'models','memristor_models','yakopic',model_syntax_ver);
        elseif strcmpi(model_ver,'picket')
            models_folder=fullfile(prj_dir,'models','memristor_models','picket',model_syntax_ver);
        elseif strcmpi(model_ver,'laiho_biolek')
            models_folder=fullfile(prj_dir,'models','memristor_models','laiho_biolek',model_syntax_ver);
        elseif strcmpi(model_ver,'univ_michigan')
            models_folder=fullfile(prj_dir,'models','memristor_models','univ_michigan',model_syntax_ver);
        elseif strcmpi(model_ver,'lineal')
            models_folder=fullfile(prj_dir,'models','memdiode','lineal',model_syntax_ver);
        elseif strcmpi(model_ver,'pcm')
            models_folder=fullfile(prj_dir,'models','PCM','1era_version',model_syntax_ver);    
        end
        
        if strcmpi(in_polarity,'unipolar')
            transfer_fcn='logsig';
        elseif strcmpi(in_polarity,'bipolar')
            transfer_fcn='tansig';
        end
        
        if ~isempty(strfind(memdiode_model,'_var'))
            fid = fopen(fullfile(models_folder,strcat(memdiode_model,'.sp')));
            tline = fgetl(fid);
            while ischar(tline)
                %disp(tline)
                if strfind(tline, '.param beta=')
                    beta=str2double(extractAfter(tline,'='));
                end
                if ~isempty(strfind(tline, '.param imax=')) || ~isempty(strfind(tline, '+ imax='))
                    Imax=str2double(extractAfter(tline,'='));
                end
                if ~isempty(strfind(tline, '.param imin=')) || ~isempty(strfind(tline, '+ imin='))
                    Imin=str2double(extractAfter(tline,'='));
                end        
                if ~isempty(strfind(tline, '.param alphamin=')) || ~isempty(strfind(tline, '+ alphamin='))
                    alphamin=str2double(extractAfter(tline,'='));
                end  
                if ~isempty(strfind(tline, '.param alphamax=')) || ~isempty(strfind(tline, '+ alphamax='))
                    alphamax=str2double(extractAfter(tline,'='));
                end        
                tline = fgetl(fid);                        
            end
            fclose(fid);
        end
    
    %% Netlist generation

        %status_guiObject.p=uipanel;
        if guiStatus_data
            %status_guiObject.log=uicontrol('Style','listbox','Enable','off','String',str_to_display,'FontName','Monospaced','Position',status_guiObject.log.Position);

            str_to_display{end+1}=sprintf("%sCreating circuit netlist for a Neural Network of %d by %d...",str_header,neu_per_layer(1), neu_per_layer(2));
            figure(status_gui);
            status_guiObject.log.String=str_to_display;
            status_guiObject.log.Value=length(str_to_display)+1;
            drawnow
        end
        
        fprintf('---> Creating circuit netlist for a Neural Network of ');
        for layer_str_i=1:length(neural_layers)-1
            fprintf('%d by %d ',neu_per_layer(layer_str_i),neu_per_layer(layer_str_i+1));
            fprintf('by ');
        end
        fprintf('...');        

        netlist_name='MLP_';
        for layer_str_i=1:length(neural_layers)-1
            netlist_name=strcat(netlist_name,num2str(neu_per_layer(layer_str_i)),'by',num2str(neu_per_layer(layer_str_i+1)),'_');
        end
        netlist_name=strcat(netlist_name,'closed_loop.net');    

        tStart_netlist=tic;
        fileID = fopen(fullfile(sim_folder,netlist_name),'w');
        fprintf(fileID,'*****************************************************************\n');
        fprintf(fileID,'*****************************************************************\n');
        fprintf(fileID,'*****************************************************************\n');
        fprintf(fileID,'*****************************************************************\n');
        fprintf(fileID,'\n');
        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
            fprintf(fileID,'.global gnd! vdd! vss!\n');
            fprintf(fileID,'.option search=''../../../models/std_cells/ibm13rfrtv/models/''\n');                                   
            %fprintf(fileID,'.option PARHIER = LOCAL\n');
            %fprintf(fileID,'.option PORT_VOLTAGE_SCALE_TO_2X = 1\n');
            %fprintf(fileID,'.option ARTIST=0 PSF=2\n');
            %fprintf(fileID,'.option PROBE=1\n');
            fprintf(fileID,'.option INGOLD=1\n');
            fprintf(fileID,'.option RUNLVL=%d\n',RUNLVL);
            fprintf(fileID,'.option OPFILE=1\n');
            fprintf(fileID,'.option SPLIT_DP=1\n');
            if strcmp(analysis,'transient')
                fprintf(fileID,'.option MEASFORM=1\n');
            else
             	fprintf(fileID,'.option MEASFORM=1\n');
            end
            if strcmpi(use_caps,'yes') && strcmpi(simEngine,'FineSim')
                fprintf(fileID,'.option finesim_fcapmin=1e-20\n');    
                fprintf(fileID,'.option finesim_resmax=1e15\n');    
                fprintf(fileID,'.option finesim_spred=0\n');
                fprintf(fileID,'.option finesim_double_precision_output=1\n');
                fprintf(fileID,'.option measdgt=8\n');
                fprintf(fileID,'.option finesim_delmax=1e-6\n');
            end
            if strcmpi(noise_transient,'yes')
                fprintf(fileID,'.option finesim_mode=spicemd\n');
            end
            fprintf(fileID,'.option CSDF=1\n');    
            fprintf(fileID,'.option SYMB=1\n');
            fprintf(fileID,'.option RELTOL=0.01\n');
            fprintf(fileID,'.option METHOD=TRAP\n');
            fprintf(fileID,'.option DELMAX=10u\n');
            fprintf(fileID,'.option LIST=1\n');
            fprintf(fileID,'.option LIS_NEW=1\n'); 
%           fprintf(fileID,'.option INTERP=1\n'); 
            %fprintf(fileID,'.option MEASFILE=1\n');
            %fprintf(fileID,'.option ITRPRT=1\n');
            fprintf(fileID,'.option POST=CSDF\n');
            if strcmp(simEngine,'FineSim')
                if strcmp(output_fmt,'binary')
                    fprintf(fileID,'.option finesim_print_to_probe=1\n');
                    fprintf(fileID,'.option finesim_output=fsdb\n');
                    fprintf(fileID,'.option finesim_vprbtol=%.12e\n', print_V_tolerance);
                    fprintf(fileID,'.option finesim_iprbtol=%.12e\n', print_I_tolerance);
                    %fprintf(fileID,'.option finesim_convlevel=1\n');      %cambie esto comentandolo 07/02/21: Tiene que estar comentado, se ve que agrega alguna conductancia en paralelo y compromete la convergencia del modelo nuevo
                    %fprintf(fileID,'.option finesim_delmax=1e-6\n');
                    %fprintf(fileID,'.option finesim_write_instance_table=1\n');
                    %fprintf(fileID,'.option finesim_write_mcparam=1\n');
                end
                if Rs<0.1
                    fprintf(fileID,'.option finesim_resmin==%.3e\n',min([Rs/2,R_in_TIA/2,Rs_PAD/2,R_sense/2]));   
                end
            elseif strcmp(simEngine,'HSPICE')
                if strcmp(output_fmt,'binary')
                    fprintf(fileID,'.option WDF=1\n');
                end
            end
            fprintf(fileID,'.temp 25\n');
            fprintf(fileID,'.lib ''allModels.inc'' tt\n');
            fprintf(fileID,'.include ''design.inc''\n');  
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/ibm13rfrvt_XOR2X4TF.sp"\n'));
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/ibm13rfrvt_DFFSRX4TF.sp"\n'));
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/ibm13rfrvt_BUFX4TF.sp"\n')); 
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/ibm13rfrvt_INVX4TF.sp"\n'));             
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/ibm13rfrvt_OR2X4TF.sp"\n'));
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/ibm13rfrvt_AND2X4TF.sp"\n'));             
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/ibm13rfrvt_AND3X4TF.sp"\n'));             
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/ibm13rfrvt_AND4X4TF.sp"\n'));             
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/ibm13rfrvt_neuron.sp"\n'));    
            fprintf(fileID,strcat('.inc "',root_str,'/models/std_cells/ibm13rfrtv/passgate_C.sp"\n'));
        end        
                
        if ispc
            fprintf(fileID,strcat('.inc "',root_str,'\\models\\memdiodo2.net"\n'));
        elseif isunix
            if strcmpi(model_ver,'DMM')
                fprintf(fileID,strcat('.inc "',root_str,filesep,'models',filesep,'memdiode',filesep,'DMM',filesep,model_syntax_ver,filesep,memdiode_model,'.sp"\n'));                
            elseif strcmpi(model_ver,'QMM')
                fprintf(fileID,strcat('.inc "',root_str,filesep,'models',filesep,'memdiode',filesep,'QMM',filesep,model_syntax_ver,filesep,memdiode_model,'.sp"\n'));                                    
            elseif strcmpi(model_ver,'QMM_SBSF')
                fprintf(fileID,strcat('.inc "',root_str,filesep,'models',filesep,'memdiode',filesep,'QMM_SBSF',filesep,model_syntax_ver,filesep,memdiode_model,'.sp"\n'));                                    
            elseif strcmpi(model_ver,'yakopic')
                fprintf(fileID,strcat('.inc "',root_str,filesep,'models',filesep,'memristor_models',filesep,'yakopic',filesep,model_syntax_ver,filesep,memdiode_model,'.sp"\n'));                                    
            elseif strcmpi(model_ver,'picket')
                fprintf(fileID,strcat('.inc "',root_str,filesep,'models',filesep,'memristor_models',filesep,'picket',filesep,model_syntax_ver,filesep,memdiode_model,'.sp"\n'));                                    
            elseif strcmpi(model_ver,'laiho_biolek')
                fprintf(fileID,strcat('.inc "',root_str,filesep,'models',filesep,'memristor_models',filesep,'laiho_biolek',filesep,model_syntax_ver,filesep,memdiode_model,'.sp"\n'));                                    
            elseif strcmpi(model_ver,'univ_michigan')
                fprintf(fileID,strcat('.inc "',root_str,filesep,'models',filesep,'memristor_models',filesep,'univ_michigan',filesep,model_syntax_ver,filesep,memdiode_model,'.sp"\n'));                                    
            elseif strcmpi(model_ver,'lineal')
                fprintf(fileID,strcat('.inc "',root_str,filesep,'models',filesep,'memdiode',filesep,'lineal',filesep,model_syntax_ver,filesep,memdiode_model,'.sp"\n'));
            elseif strcmpi(model_ver,'PCM')
                fprintf(fileID,strcat('.inc "',root_str,filesep,'models',filesep,'PCM',filesep,'1era_version',filesep,model_syntax_ver,filesep,memdiode_model,'.sp"\n'));                
            end
        end
        fprintf(fileID,'\n'); 
        
        
%% This part of the code creates the subcircuit of the crossabr structures
        memristor_param_lst = rram_crossbar(fileID, simEngine, neural_layers, neu_per_layer_sub, connections, memdiode_model, use_caps, dev_polarity);

%% This parte of the code creates the subcircuit of a selector     
        selector(fileID,simEngine,sel_type,R_off_value)

%% This parte of the code creates the subcircuit of a non-overlapping clk generator        
        nonOv_clk(fileID,simEngine)

%% This parte of the code creates an ideal TIA of gain 1000     
        TIA_circuit(fileID,simEngine)
        
%% This parte of the code creates the logic gates necesary for the circuit 

        and1x4tf(fileID,simEngine) 
        and5x4tf(fileID,simEngine)         
        and6x4tf(fileID,simEngine) 
        and7x4tf(fileID,simEngine) 
        and8x4tf(fileID,simEngine)    
        and9x4tf(fileID,simEngine) 
        and10x4tf(fileID,simEngine) 
        and12x4tf(fileID,simEngine) 
        and15x4tf(fileID,simEngine)    
        and16x4tf(fileID,simEngine)    
        and18x4tf(fileID,simEngine)     
        and28x4tf(fileID,simEngine) 
        and32x4tf(fileID,simEngine) 
        demux4tf(fileID,simEngine)
        
%% This parte of the code creates the subcircuit of a N-bits decoder.
        nBit_decoder(fileID,simEngine,neu_per_layer_sub)

%% This parte of the code creates the subcircuit of a N-bits counter.   
        nBit_counter(fileID,simEngine,neu_per_layer_sub,max_pulses)
        
%% This parte of the code creates the neuron subcircuit 
        lin_neuron(fileID, ...
                   partitions_N, ...
                   'rSense',R_sense, ...
                   'tiaGain',TIA_gain,...
                   'TIA',TIA, ...
                   'Gsense',Gsense,...
                   'simEngine',simEngine);

%% This parte of the code creates the neuron subcircuit for the logsig neuron
        log_neuron(fileID, ...
                   partitions_N, ...
                   'rSense',R_sense, ...
                   'tiaGain',TIA_gain,...
                   'TIA',TIA, ...
                   'Gsense',Gsense,...
                   'simEngine',simEngine, ...
                   'neurons',neurons, ...
                   'readVoltage', read_voltage);    

%% This parte of the code creates the neuron subcircuit for the tansig neuron
         tan_neuron(fileID, ...
                   partitions_N, ...
                   'rSense',R_sense, ...
                   'tiaGain',TIA_gain,...
                   'TIA',TIA, ...
                   'Gsense',Gsense,...
                   'simEngine',simEngine, ...
                   'neurons', neurons, ...
                   'readVoltage', read_voltage);               
        
%% This parte of the code creates the subcircuit of a N-bits decoder.            
        fprintf(fileID,'*#########################   %d layer synchronizer SUBCKT  ############################\n',length(neural_layers)-1);
        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
            fprintf(fileID,'.subckt sync_W_OK_%d_layers ', length(neural_layers)-1);
            for i=1:length(neural_layers)-1
                fprintf(fileID,'WOK_layer_%d ', i);    
            end
            fprintf(fileID,'output clk PRE_1 clear\n\n');
            for i=1:length(neural_layers)-1
                fprintf(fileID,'Xdff%d clk_en gnd! output_%d RQn%d clearn WOK_layer_%dn dffsrx4tf\n',i , i, i, i);
                fprintf(fileID,'Xinv%d WOK_layer_%d WOK_layer_%dn invx4tf\n',i, i ,i);
                fprintf(fileID,'RQn%d RQn%d gnd! R=1Meg\n\n',i , i);
            end
            fprintf(fileID,'Xinv_clear clear clearn invx4tf\n\n');
            fprintf(fileID,'Xand%d ',length(neural_layers)-1);
            for i=1:length(neural_layers)-1
                fprintf(fileID,'output_%d ',i);
            end
            fprintf(fileID,'output and%dx4tf\n', length(neural_layers)-1);
            fprintf(fileID,'Xdff_gral clk output clk_en clk_en_n clearn vdd! dffsrx4tf\n');
        else

        end
        fprintf(fileID,'*########################################################################\n');
        fprintf(fileID,'.ends sync_W_OK_%d_layers\n\n', length(neural_layers)-1);
        
%% This parte of the code creates the top circuit level             
        fprintf(fileID,'*========================================================================================================================================================\n');
        fprintf(fileID,'*##############################################################  TOP CIRCUIT  ###########################################################################\n\n');
      
        power_meas_statement_idx=1;
        for CPA_i=1:length(neural_layers)-1
            W_matrix_init_CPA=W_matrix_init{1,CPA_i};
            W_matrix_CPA=W_matrix{1,CPA_i};
            fprintf(fileID,'*======================================================================================================================================\n');
            fprintf(fileID,'*Layers %s-%s, including interconections (CPA), input and output neurons and programming electronic', neural_layers{CPA_i},neural_layers{CPA_i+1});
            fprintf(fileID,'****************************************************************************************************************************************\n');  
            fprintf(fileID,'*#######################################################################################################################################\n');  
            fprintf(fileID,'\n');  
            for pos_neg_i=1:2    
                for partitions_i=1:partitions_N(CPA_i)
                    for partitions_ii=1:partitions_N(CPA_i+1)
                        fprintf(fileID,'*========================================================================\n');
                        fprintf(fileID,'*Layers (%s-%s) interconnections (CPA): Partition %d-%d, %s polarity\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                        fprintf(fileID,'*************************************************************************\n');  
        
                        fprintf(fileID,'XMEMD_network_%s-%s_P%d-%d_%s ', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                        for i=1:neu_per_layer_sub(CPA_i)
                            fprintf(fileID,'%s%s%d_PAD_1_P%d-%d_%s ', neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                        end
                        for i=1:neu_per_layer_sub(CPA_i+1)
                            fprintf(fileID,'%s%s%d_PAD_1_P%d-%d_%s ', neural_layers{CPA_i+1}, neural_layers{CPA_i}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i});    
                        end
                        if isempty(strfind(memdiode_model,'_var'))
                            fprintf(fileID,'MEMD_network_%d_by_%d ', neu_per_layer_sub(CPA_i), neu_per_layer_sub(CPA_i+1));
                            for i=1:neu_per_layer_sub(CPA_i)
                                for ii=1:neu_per_layer_sub(CPA_i+1)
                                    fprintf(fileID,'H0_r%d_c%d=%f ', i, ii, W_matrix_init_CPA(i+(partitions_i-1)*neu_per_layer_sub(CPA_i),ii+(partitions_ii-1)*neu_per_layer_sub(CPA_i+1),pos_neg_i));    
                                end
                            end                    
                        else
                            fprintf(fileID,'MEMD_network_%d_by_%d \n', neu_per_layer_sub(CPA_i), neu_per_layer_sub(CPA_i+1));
                            for i=1:neu_per_layer_sub(CPA_i)
                                for ii=1:neu_per_layer_sub(CPA_i+1)
                                    fprintf(fileID,'+ H0_r%d_c%d=%f Imax_A%d_B%d=%.4e Imin_A%d_B%d=%.4e \n',i, ii, W_matrix_init_CPA(i+(partitions_i-1)*neu_per_layer_sub(1),ii,pos_neg_i),i, ii,Imax*normrnd(1,sigmaImax),i, ii,Imin*normrnd(1,sigmaImin));    
                                end
                            end
                        end
                        fprintf(fileID,'\n');

                        upper_layer_RS_lst=memristor_param_lst{CPA_i,1}.upper_layer_RS_lst;
                        lower_layer_RS_lst=memristor_param_lst{CPA_i,1}.lower_layer_RS_lst;
                        memdiode_Rcs_lst=memristor_param_lst{CPA_i,1}.memdiode_Rcs_lst;
                        memdiode_RS_lst=memristor_param_lst{CPA_i,1}.memdiode_RS_lst;

                        curr_partition_str=sprintf('xmemd_network_%s-%s_P%d-%d_%s', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                        for upper_layer_RS_lst_index=1:length(upper_layer_RS_lst)
                            power_meas_statement{power_meas_statement_idx,1}=sprintf('.meas TRAN %s RMS p(%s.%s)',sprintf('CPA_%s%d-%d_%s-%s.RU%d',pos_neg{pos_neg_i},partitions_i,partitions_ii,neural_layers{CPA_i}, neural_layers{CPA_i+1},upper_layer_RS_lst_index), curr_partition_str, upper_layer_RS_lst{upper_layer_RS_lst_index,1});
                            power_meas_statement_idx=power_meas_statement_idx+1;
                        end
                        for lower_layer_RS_lst_index=1:length(lower_layer_RS_lst)
                            power_meas_statement{power_meas_statement_idx,1}=sprintf('.meas TRAN %s RMS p(%s.%s)',sprintf('CPA_%s%d-%d_%s-%s.RL%d',pos_neg{pos_neg_i},partitions_i,partitions_ii,neural_layers{CPA_i}, neural_layers{CPA_i+1},lower_layer_RS_lst_index), curr_partition_str, lower_layer_RS_lst{lower_layer_RS_lst_index,1});
                            power_meas_statement_idx=power_meas_statement_idx+1;
                        end
                        for memdiode_Rcs_lst_index=1:length(memdiode_Rcs_lst)
                            power_meas_statement{power_meas_statement_idx,1}=sprintf('.meas TRAN %s RMS p(%s.%s)',sprintf('CPA_%s%d-%d_%s-%s.RC%d',pos_neg{pos_neg_i},partitions_i,partitions_ii,neural_layers{CPA_i}, neural_layers{CPA_i+1},memdiode_Rcs_lst_index), curr_partition_str, memdiode_Rcs_lst{memdiode_Rcs_lst_index,1});
                            power_meas_statement_idx=power_meas_statement_idx+1;
                        end
                        for memdiode_RS_lst_index=1:length(memdiode_RS_lst)
                            power_meas_statement{power_meas_statement_idx,1}=sprintf('.meas TRAN %s RMS PAR(''(V(%s.%s)-V(%s.%s))*I(%s.%s)'')',sprintf('CPA_%s%d-%d_%s-%s.MD%d',pos_neg{pos_neg_i},partitions_i,partitions_ii, neural_layers{CPA_i}, neural_layers{CPA_i+1}, memdiode_RS_lst_index), curr_partition_str, memdiode_RS_lst{memdiode_RS_lst_index,1}, curr_partition_str, memdiode_RS_lst{memdiode_RS_lst_index,2}, curr_partition_str, memdiode_RS_lst{memdiode_RS_lst_index,3});
                            power_meas_statement_idx=power_meas_statement_idx+1;
                        end
                        
                         %% Drivers of the first layer
                        fprintf(fileID,'\n');
                        fprintf(fileID,'*========================================================================\n');
                        fprintf(fileID,'*Drivers for the input neurons of intercon. layer %s-%s (Neural layer %s, CPA: %d-%d-%s)\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                        fprintf(fileID,'*************************************************************************\n');
                        for i=1:neu_per_layer_sub(CPA_i)
                            time_accu=0;
                            fprintf(fileID,'Xselector_RW_%s%s%d-1-P%d-%d_%s %s%s%d_PAD_1_P%d-%d_%s RW_%s-%s write_signal_%s%s%d_1_P%d-%d_%s Neuron%s%d_P%d-1 selector\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, i, partitions_i);    
                            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                fprintf(fileID,'Xselector_RW_%s%s%d-2-P%d-%d_%s write_signal_%s%s%d_1_P%d-%d_%s Vcont%s%s%d-1_out_P%d-%d_%s Vpulsedsignal RHZ_%s%s%d_P%d-%d_%s selector\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                                %fprintf(fileID,'RHZ_%s%d RHZ_%s%d gnd! R=''HZ''\n',neural_layers{1}, i,neural_layers{1}, i);
                                fprintf(fileID,'VHZ_%s%s%d_P%d-%d_%s RHZ_%s%s%d_P%d-%d_%s gnd! dc=%f\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, unselected_row_voltage);                                  
                                fprintf(fileID,'Xand%s%s%d_P%d-%d_%s Vcont%s-%s%d-1 enable_write_%s-%s_P%d-%d_%s Vcont%s%s%d-1_out_P%d-%d_%s and2x4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i});                                
                                if strcmpi(pos_neg{pos_neg_i},'pos') && partitions_ii==1 && CPA_i==1

                                    fprintf(fileID,'VNeuron%s%d_P%d Neuron%s%d_P%d-1 gnd! PWL PWLFILE=''%s/input_sources/V_input_%s%d_P%d.csv''\n', neural_layers{CPA_i}, i, partitions_i, neural_layers{CPA_i}, i, partitions_i, sim_folder, neural_layers{CPA_i}, i, partitions_i);
                                    if ~exist(fullfile(sim_folder,'input_sources'),'dir')
                                        mkdir(fullfile(sim_folder,'input_sources'))
                                    end
                                    fileID_vsources=fopen(fullfile(sim_folder,'input_sources',sprintf('V_input_%s%d_P%d.csv', neural_layers{CPA_i}, i, partitions_i)),'w');
                                    fprintf(fileID_vsources,'0 0 %.12e 0\n',read_delay);
                                    time_accu=time_accu+read_delay;
                                    for nVectors=1:size(inputMatrix,2)
                                        fst_time=time_accu+in_vect_r_f_time;
                                        last_time=fst_time+1/input_vector_freq;
                                        fprintf(fileID_vsources,'%.12e %f %.12e %f\n',fst_time,inputMatrix(i+(partitions_i-1)*neu_per_layer_sub(1),nVectors)*read_voltage,last_time,inputMatrix(i+(partitions_i-1)*neu_per_layer_sub(1),nVectors)*read_voltage);
                                        time_accu=last_time;
                                    end
                                    time_accu=time_accu+in_vect_r_f_time;
                                    if sim_WR==0
                                        fprintf(fileID_vsources,'%.12e %f %.12e %f',time_accu, inputMatrix(i+(partitions_i-1)*neu_per_layer_sub(1),nVectors)*read_voltage, simulation_time, inputMatrix(i+(partitions_i-1)*neu_per_layer_sub(1),nVectors)*read_voltage);
                                    else
                                        fprintf(fileID_vsources,'%.12e %f %.12e %f',time_accu, inputMatrix(i+(partitions_i-1)*neu_per_layer_sub(1),nVectors)*read_voltage, time_accu*1.05, inputMatrix(i+(partitions_i-1)*neu_per_layer_sub(1),nVectors)*read_voltage);
                                    end
                                    fclose(fileID_vsources);
                                end
                            else
                                fprintf(fileID,'Xselector_RW_%s%d-2-%d write_signal_%s%d_1_%d Vcont%s%d-1_out_%d Vpulsedsignal RHZ_%s%d selector\n', neural_layers{1}, i, pos_neg{pos_neg_i}, neural_layers{1}, i, pos_neg{pos_neg_i}, neural_layers{1}, i, pos_neg{pos_neg_i}, neural_layers{1}, i);
                                fprintf(fileID,'Aand%s%d Vcont%s%d-1 enable_write VDD VDD VDD VDD Vcont%s%d-1_out_%d 0 AND\n',neural_layers{1}, i,neural_layers{1}, i, neural_layers{1}, i);   
                                %fprintf(fileID,'RHZ_%s%d RHZ_%s%d 0 {HZ}\n',neural_layers{1}, i,neural_layers{1}, i);
                                fprintf(fileID,'VHZ_%s%d RHZ_%s%d 0 %f\n',neural_layers{1}, i,neural_layers{1}, i,0);                                  
                                if strcmpi(pos_neg{pos_neg_i},'pos')
                                    fprintf(fileID,'VNeuron%s%d Neuron%s%d-1 0 PWL(0 0 %.12e 0 ', neural_layers{1}, i, neural_layers{1}, i, read_delay);
                                    time_accu=time_accu+read_delay;
                                    for nVectors=1:size(inputMatrix,2)
                                        fst_time=time_accu+in_vect_r_f_time;
                                        last_time=fst_time+1/input_vector_freq;
                                        fprintf(fileID,'%.12e %f %.12e %f ',fst_time,inputMatrix(i,nVectors)*read_voltage,last_time,inputMatrix(i,nVectors)*read_voltage);
                                        time_accu=last_time;
                                    end
                                    time_accu=time_accu+in_vect_r_f_time;
                                    fprintf(fileID,'%.12e %f %.12e %f)\n',time_accu, read_voltage, simulation_time, read_voltage);
                                end
                            end
                            fprintf(fileID,'*************************************************************************\n');           
                        end
                        fprintf(fileID,'\n');

                        %% Drivers of the second layer
                        fprintf(fileID,'*========================================================================\n');
                        fprintf(fileID,'*Drivers for the output neurons of intercon. layer %s-%s (Neural layer %s, CPA: %d-%d-%s)\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                        fprintf(fileID,'*************************************************************************\n');
                        for i=1:neu_per_layer_sub(CPA_i+1)  
                            fprintf(fileID,'Xselector_RW_%s%s%d-1-P%d-%d_%s %s%s%d_PAD_1_P%d-%d_%s RW_%s-%s write_signal_%s%s%d_1_P%d-%d_%s Neuron%s%d_P%d-%d-%s selector\n', neural_layers{CPA_i+1}, neural_layers{CPA_i}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i});                        
                            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                fprintf(fileID,'Xselector_RW_%s%s%d-2-P%d-%d_%s write_signal_%s%s%d_1_P%d-%d_%s Vcont%s%s%d-1_out_P%d-%d_%s comparator_%s-%s_P%d-%d_%s RHZ_%s%s%d_P%d-%d_%s selector\n', neural_layers{CPA_i+1}, neural_layers{CPA_i}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                                fprintf(fileID,'VHZ_%s%s%d_P%d-%d_%s RHZ_%s%s%d_P%d-%d_%s gnd! dc=%f\n',neural_layers{CPA_i+1}, neural_layers{CPA_i}, i,  partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, i,  partitions_i, partitions_ii, pos_neg{pos_neg_i}, unselected_col_voltage);                  
                                fprintf(fileID,'Xand%s%s%d_P%d-%d_%s Vcont%s-%s%d-1 enable_write_%s-%s_P%d-%d_%s Vcont%s%s%d-1_out_P%d-%d_%s and2x4tf\n',neural_layers{CPA_i+1}, neural_layers{CPA_i}, i,  partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, i, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, i,  partitions_i, partitions_ii, pos_neg{pos_neg_i});                                
                                if strcmpi(Gsense,'yes') || CPA_i==length(neural_layers)-1
                                    fprintf(fileID,'RNeuron%s%d_P%d-%d_%s Neuron%s%d_P%d-%d-%s gnd! R=100Meg\n', neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i});           
                                else
                                    fprintf(fileID,'RNeuron%s%d_P%d-%d_%s Neuron%s%d_P%d-%d-%s gnd! R=100Meg\n', neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i});           
                                end
                            else
                                fprintf(fileID,'Sselector_RW_%s%d-2-P%d_%s write_signal_%s%d_1_P%d_%s Vcont%s%d-1_out_P%d_%s comparator_P%d_%s RHZ_%s%d_P%d_%s selector\n', neural_layers{2}, i, partitions_i, pos_neg{pos_neg_i}, neural_layers{2}, i, partitions_i, pos_neg{pos_neg_i}, neural_layers{2}, i, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, neural_layers{2}, i, partitions_i, pos_neg{pos_neg_i});
                                fprintf(fileID,'VHZ_%s%d_P%d_%s RHZ_%s%d_P%d_%s 0 dc=%f\n',neural_layers{2}, i,  partitions_i, pos_neg{pos_neg_i}, neural_layers{2}, i,  partitions_i, pos_neg{pos_neg_i}, uselected_col_voltage);                  
                                fprintf(fileID,'Aand%s%d_P%d_%s Vcont%s%d-1 enable_write_P%d_%s VDD VDD VDD VDD Vcont%s%d-1_out_P%d_%s 0 AND\n',neural_layers{2}, i,  partitions_i, pos_neg{pos_neg_i}, neural_layers{2}, i,  partitions_i, pos_neg{pos_neg_i}, neural_layers{2}, i,  partitions_i, pos_neg{pos_neg_i});                                
                                fprintf(fileID,'RNeuron%s%d_P%d_%s Neuron%s%d_P%d-%s 0 1\n', neural_layers{2}, i, partitions_i, pos_neg{pos_neg_i}, neural_layers{2}, i, partitions_i, pos_neg{pos_neg_i});      
                            end
                            fprintf(fileID,'*************************************************************************\n');    
                        end 
                        fprintf(fileID,'\n');
                        if sim_WR==1
                            fprintf(fileID,'*========================================================================\n');                        
                            fprintf(fileID,'*Voltage references for memdiode programation of layers %s-%s (CPA: %d-%d-%s)\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            fprintf(fileID,'*************************************************************************\n');
                            for i=1:neu_per_layer_sub(CPA_i)
                                for ii=1:neu_per_layer_sub(CPA_i+1)
                                    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                        fprintf(fileID,'Vref_%s%d-%s%d_P%d-%d_%s ref_%s%d-%s%d-2_P%d-%d_%s gnd! dc=%f\n', neural_layers{CPA_i}, i, neural_layers{CPA_i+1}, ii, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, i, neural_layers{CPA_i+1}, ii, partitions_i, partitions_ii, pos_neg{pos_neg_i}, W_matrix_CPA(i+(partitions_i-1)*neu_per_layer_sub(CPA_i),ii+(partitions_ii-1)*neu_per_layer_sub(CPA_i+1),pos_neg_i));        
                                        fprintf(fileID,'gi%s%d-%s%d_P%d-%d_%s ref_%s%d-%s%d-1_P%d-%d_%s ref_%s%d-%s%d-2_P%d-%d_%s vcr pwl(1) Vcont%s-%s%d-1 gnd! 0.495,%s 0.505,1\n', neural_layers{CPA_i}, i, neural_layers{CPA_i+1}, ii, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, i, neural_layers{CPA_i+1}, ii, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, i, neural_layers{CPA_i+1}, ii, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, i, R_off_value);
                                        fprintf(fileID,'gi%s%d-%s%d_P%d-%d_%s V_ref_%s-%s_P%d-%d_%s ref_%s%d-%s%d-1_P%d-%d_%s vcr pwl(1) Vcont%s-%s%d-1 gnd! 0.495,%s 0.505,1\n', neural_layers{CPA_i+1}, ii, neural_layers{CPA_i}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, i, neural_layers{CPA_i+1}, ii, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, ii, R_off_value);                
                                    else
                                        %fprintf(fileID,'Vref_%s%d-%s%d_%d ref_%s%d-%s%d-2_%d 0 %f\n', neural_layers{1}, i, neural_layers{2}, ii, pos_neg, neural_layers{1}, i, neural_layers{2}, ii, pos_neg, W_matrix(i,ii,pos_neg));                            
                                        %fprintf(fileID,'S%s%d-%s%d_%d ref_%s%d-%s%d-1_%d ref_%s%d-%s%d-2_%d Vcont%s%d-1 0 SW1\n', neural_layers{1}, i, neural_layers{2}, ii, pos_neg, neural_layers{1}, i, neural_layers{2}, ii, pos_neg, neural_layers{1}, i, neural_layers{2}, ii, pos_neg, neural_layers{1}, i);
                                        %fprintf(fileID,'S%s%d-%s%d_%d V_ref_%d ref_%s%d-%s%d-1_%d Vcont%s%d-1 0 SW1\n', neural_layers{2}, ii, neural_layers{1}, i, pos_neg, pos_neg, neural_layers{1}, i, neural_layers{2}, ii, pos_neg, neural_layers{2}, ii);
                                    end
                                end
                            end
                            fprintf(fileID,'*************************************************************************\n');   
                            fprintf(fileID,'\n');
                        end
                        %% Event detectors for each crossbar array
                        fprintf(fileID,'*========================================================================\n');
                        fprintf(fileID,'*Events (write) detection circuit for CPA %s-%s partition array %d-%d, group %s \n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                        fprintf(fileID,'*************************************************************************\n');
                        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim') %HSPICE piece of code
                            if sim_WR==1
                                %Esto es usando TIA
                                %fprintf(fileID,'Xtia_%s-%s_P%d-%d_%s comparator_%s-%s_P%d-%d_%s comparator_%s-%s_P%d-%d_%s-1 TIA_ideal\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                                fprintf(fileID,'E1_%s-%s_P%d-%d_%s event_in_%s-%s_P%d-%d_%s! gnd! vol=''((V(comparator_%s-%s_P%d-%d_%s-1)-V(V_ref_%s-%s_P%d-%d_%s))*2e6)'' min=0 max=1.2\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});

                                %Esto es sensando tensión
                                %fprintf(fileID,'E1_%s-%s_P%d-%d_%s event_in_%s-%s_P%d-%d_%s! gnd! vol=''((V(comparator_%s-%s_P%d-%d_%s-1)-V(V_ref_%s-%s_P%d-%d_%s))*2e6)'' min=0 max=1.2\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            else
                                fprintf(fileID,'E1_%s-%s_P%d-%d_%s event_in_%s-%s_P%d-%d_%s! gnd! vol=''((V(comparator_%s-%s_P%d-%d_%s-1)-0)*2e6)'' min=0 max=1.2\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            end
                            fprintf(fileID,'gi_%s-%s_P%d-%d_%s comparator_%s-%s_P%d-%d_%s comparator_%s-%s_P%d-%d_%s-1 vcr pwl(1) enable_read_%s-%s_P%d-%d_%s gnd! 0.1,%s 0.9,1\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, R_off_value);
                            
                            fprintf(fileID,'Xsync_read_enable_%s-%s_P%d-%d_%s sys_clk enable_write_%s-%s_P%d-%d_%s control_enable_read_%s-%s_P%d-%d_%s control_enable_read_n_%s-%s_P%d-%d_%s sr_rw_sync vdd! dffsrx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            fprintf(fileID,'Xand_read_enable_%s-%s_P%d-%d_%s control_enable_read_%s-%s_P%d-%d_%s enable_read_pulse enable_read_%s-%s_P%d-%d_%s and2x4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            
                            fprintf(fileID,'Raux1_%s-%s_P%d-%d_%s comparator_%s-%s_P%d-%d_%s-1 gnd! R=1Meg\n',neural_layers{CPA_i+1}, neural_layers{CPA_i}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            %fprintf(fileID,'Raux2_%s-%s_P%d-%d_%s comparator_%s-%s_P%d-%d_%s-1 gnd! R=1Meg\n',neural_layers{CPA_i+1}, neural_layers{CPA_i}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            fprintf(fileID,'Rshunt_%s-%s_P%d-%d_%s comparator_%s-%s_P%d-%d_%s gnd! R=''R_shunt''\n',neural_layers{CPA_i+1}, neural_layers{CPA_i}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            %fprintf(fileID,'Cshunt_%s-%s_P%d-%d_%s comparator_%s-%s_P%d-%d_%s gnd! C=1p\n',neural_layers{CPA_i+1}, neural_layers{CPA_i}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            fprintf(fileID,'\n');
                            fprintf(fileID,'Xbuf_event_%s-%s_P%d-%d_%s event_in_%s-%s_P%d-%d_%s! event_in_m_%s-%s_P%d-%d_%s! bufx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            fprintf(fileID,'Xinv_event_%s-%s_P%d-%d_%s event_in_m_%s-%s_P%d-%d_%s! event_out_%s-%s_P%d-%d_%s! invx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii,  pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});   
                            if strcmpi(pos_neg{pos_neg_i},'pos') && partitions_i==1 && partitions_ii==1 && CPA_i==1
                                fprintf(fileID,'XinvReset resetDFFs sr_rw_sync invx4tf\n');            
                            end
                            fprintf(fileID,'Xand_event_TO__%s-%s_P%d-%d_%s time_out_n_%s-%s event_out_%s-%s_P%d-%d_%s! event_out_TO_%s-%s_P%d-%d_%s! and2x4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            fprintf(fileID,'Xevent_det_%s-%s_P%d-%d_%s clk_event_%s-%s gnd! EVENT_%s-%s_P%d-%d_%s! EVENTn_%s-%s_P%d-%d_%s! sr_rw_sync event_out_TO_%s-%s_P%d-%d_%s! dffsrx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});                 %D=gnd! clk=sys_clk Q=EVENT! Qn=EVENTn! rs=clear_n ss=event_out!
                            fprintf(fileID,'Xenable_write_%s-%s_P%d-%d_%s sys_clk EVENTn_%s-%s_P%d-%d_%s! enable_write_%s-%s_P%d-%d_%s enable_write_n_%s-%s_P%d-%d_%s EVENTn_%s-%s_P%d-%d_%s! vdd! dffsrx4tf\n',neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});            
                        else                         %LTspice piece of code
                            fprintf(fileID,'XU1_P%d_%s comparator_P%d_%s-1 V_ref_P%d_%s VDD VSS event_in_P%d_%s! level.2 Avol=1000Meg GBW=1000Meg Slew=1000Meg ilimit=25m rail=0 Vos=0 phimargin=45 en=0 enk=0 in=0 ink=0 Rin=500Meg\n', partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i});
                            fprintf(fileID,'Si_P%d_%s comparator_P%d_%s comparator_P%d_%s-1 enable_read_pulse 0 SW1\n',  partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i});
                            fprintf(fileID,'Raux_%s_P%d_%s comparator_P%d_%s-1 gnd! 10k\n',neural_layers{2}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i});
                            fprintf(fileID,'Rshunt_%s_P%d_%s comparator_P%d_%s gnd! {R_shunt}\n',neural_layers{2}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i});
                            fprintf(fileID,'\n');
                            fprintf(fileID,'Abuf_event_P%d_%s event_in_P%d_%s! 0 0 0 0 0 event_in_m_P%d_%s! 0 BUF\n', partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i});
                            fprintf(fileID,'Ainv_event_P%d_%s event_in_m_P%d_%s! 0 0 0 0 0 event_out_P%d_%s! INV\n', partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i});   
                            if strcmpi(pos_neg{pos_neg_i},'pos') && partitions_i==1
                                fprintf(fileID,'XinvReset resetDFFs sr_rw_sync invx4tf\n');            
                            end
                            fprintf(fileID,'Vlow_%d V_low_%d 0 0\n', pos_neg{pos_neg_i}, pos_neg{pos_neg_i});
                            fprintf(fileID,'Aevent_det_P%d_%s V_low_%d 0 clk_event_%s-%s event_out_P%d_%s! sr_rw_sync EVENTn_P%d_%s! EVENT_P%d_%s! 0 DFLOP Vhigh=1 Vlow=0 Trise=5n Tfall=5n Td=5n\n', partitions_i, pos_neg{pos_neg_i}, partitions_i,  neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i});                 %D=gnd! clk=sys_clk Q=EVENT! Qn=EVENTn! rs=clear_n ss=event_out!
                            fprintf(fileID,'Abuf_enable_write_P%d_%s EVENTn_P%d_%s! enable_write_P%d_%s BUF\n', partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i});
                        end
                        fprintf(fileID,'*************************************************************************\n');                        
                    end
                end
            end


            %% Hidden layer of neuros
            fprintf(fileID,'*========================================================================\n');
            fprintf(fileID,'*Hidden layer of neurons (layer %s) \n', neural_layers{CPA_i+1});
            fprintf(fileID,'*************************************************************************\n');            
            if CPA_i+1<=length(neural_layers)
                for partitions_ii=1:partitions_N(CPA_i+1)
                    for i=1:neu_per_layer_sub(CPA_i+1)
                        fprintf(fileID,'XNeuron%s%d_P%d-1 ', neural_layers{CPA_i+1}, i, partitions_ii);
                        for pos_neg_i=1:2
                            for partitions_i=1:partitions_N(CPA_i)
                                fprintf(fileID,'Neuron%s%d_P%d-%d-%s ',neural_layers{CPA_i+1}, i, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                            end
                        end
                        if CPA_i+1<length(neural_layers)
                            if strcmpi(transfer_fcn,'logsig')
                                fprintf(fileID,'Neuron%s%d_P%d-1 gnd! logsig_neuron_%dx2_in\n', neural_layers{CPA_i+1}, i, partitions_ii, partitions_N(CPA_i));
                            else
                                fprintf(fileID,'Neuron%s%d_P%d-1 gnd! tansig_neuron_%dx2_in\n', neural_layers{CPA_i+1}, i, partitions_ii, partitions_N(CPA_i));
                            end
                        else
                            fprintf(fileID,'Neuron%s%d_P%d-1 gnd! neuron_%dx2_in\n', neural_layers{CPA_i+1}, i, partitions_ii, partitions_N(CPA_i));
                            fprintf(fileID,'RNeuron%s%d_P%d-1 Neuron%s%d_P%d-1 gnd! R=100Meg\n', neural_layers{CPA_i+1}, i, partitions_ii, neural_layers{CPA_i+1}, i, partitions_ii);
                        end
                    end
                end
                fprintf(fileID,'\n');
            end
            fprintf(fileID,'*************************************************************************\n');                        
           
            
            %% Counter of written rows        
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                
                fprintf(fileID,'Xdecoder%s-%s ', neural_layers{CPA_i},neural_layers{CPA_i+1});
                for i=ceil(log2(neu_per_layer_sub(CPA_i))):-1:1
                    fprintf(fileID,'S_%s-%s%d-1 ', neural_layers{CPA_i}, neural_layers{CPA_i+1}, i);
                end
                for i=1:neu_per_layer_sub(CPA_i)
                    fprintf(fileID,'Vcont%s-%s%d-1 ', neural_layers{CPA_i}, neural_layers{CPA_i+1}, i);
                end
                fprintf(fileID,'decoder_%d_output\n', neu_per_layer_sub(CPA_i));
                
                fprintf(fileID,'Xcounter_%s-%s eventsA_%s-%s write_ok_pulse!_%s-%s reset_count_deco_%s-%s ', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});
                for i=1:ceil(log2(neu_per_layer_sub(CPA_i)))
                    fprintf(fileID,'S_%s-%s%d-1 ', neural_layers{CPA_i}, neural_layers{CPA_i+1}, i);
                end                
                fprintf(fileID,'counter_%d_pulses\n',  neu_per_layer_sub(CPA_i));            
                fprintf(fileID,'Xinv_wr_pulse_%s-%s write_ok_pulse!_%s-%s write_ok_pulse_n!_%s-%s invx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});
                fprintf(fileID,'X_WR_OK_%s%s sys_clk gnd! write_ok!_%s-%s write_okn!_%s-%s resetDFFsn write_ok_pulse_n!_%s-%s dffsrx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});  
            else
                %fprintf(fileID,'Xdecoder%s ', neural_layers{1});
                %for i=1:neu_per_layer_sub(1)
                    %fprintf(fileID,'Vcont%s%d-1 ', neural_layers{1}, i);
                %end
                %fprintf(fileID,'eventsA WR_user reset_count_deco decoder_%d_output\n', neu_per_layer_sub(1));
                %fprintf(fileID,'Xcounter%s eventsA write_ok_pulse! reset_count_deco counter_%d_pulses\n', neural_layers{1}, neu_per_layer_sub(1));                                                                       %D=write_ok! clk=sys_clk Q=EVENT! Qn=EVENTn! rs=clear_n ss=event_out!                       
            end

            %% Counter of written columns
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                fprintf(fileID,'Xdecoder%s-%s ', neural_layers{CPA_i+1},neural_layers{CPA_i});
                for i=ceil(log2(neu_per_layer_sub(CPA_i+1))):-1:1
                    fprintf(fileID,'S_%s-%s%d-1 ', neural_layers{CPA_i+1}, neural_layers{CPA_i}, i);
                end                
                for i=1:neu_per_layer_sub(CPA_i+1)
                    fprintf(fileID,'Vcont%s-%s%d-1 ', neural_layers{CPA_i+1}, neural_layers{CPA_i}, i);
                end
                fprintf(fileID,'decoder_%d_output\n', neu_per_layer_sub(CPA_i+1));
                fprintf(fileID,'Xcounter_%s-%s clk_event_%s-%s eventsA_%s-%s reset_count_deco_%s-%s ', neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});
                for i=1:ceil(log2(neu_per_layer_sub(CPA_i+1)))
                    fprintf(fileID,'S_%s-%s%d-1 ', neural_layers{CPA_i+1}, neural_layers{CPA_i}, i);
                end                
                fprintf(fileID,'counter_%d_pulses\n', neu_per_layer_sub(CPA_i+1));
                fprintf(fileID,'\n');
            else
                %fprintf(fileID,'Xdecoder%s ', neural_layers{2});
                %for i=1:neu_per_layer_sub(2)
                    %fprintf(fileID,'Vcont%s%d-1 ', neural_layers{2}, i);
                %end
                %fprintf(fileID,'EVENT! WR_user reset_count_deco decoder_%d_output\n', neu_per_layer_sub(2));
                %fprintf(fileID,'Xcounter%s EVENT! eventsA reset_count_deco counter_%d_pulses\n', neural_layers{2}, neu_per_layer_sub(2));
                %fprintf(fileID,'\n');
            end
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                fprintf(fileID,'Xor%s%s%d resetDFFs RWn_%s-%s reset_count_deco_%s-%s or2x4tf\n', neural_layers{CPA_i},neural_layers{CPA_i+1}, i, neural_layers{CPA_i},neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});   
            else
    %             fprintf(fileID,'Aor%s%d resetDFFs RWn 0 0 0 0 reset_count_deco 0 OR\n',neural_layers{1}, i);   
            end
            fprintf(fileID,'\n');

            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim') %HSPICE piece of code
                for pos_neg_i=1:2
                    fprintf(fileID,'Xevent_sync_%s-%s_%s ', neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i});    
                    for partitions_i=1:partitions_N(CPA_i)
                        for partitions_ii=1:partitions_N(CPA_i+1)
                            fprintf(fileID,'event_%s-%s_P%d-%d_%s! ',neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});    
                        end
                    end
                    fprintf(fileID,'EVENT_%s-%s_%s! and%dx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i}, partitions_N(CPA_i)*partitions_N(CPA_i+1));
                end
                fprintf(fileID,'Xclk_sync_%s-%s EVENT_%s-%s_pos! EVENT_%s-%s_neg! clk_event_n_%s-%s EVENT_S_%s-%s! and3x4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});    
                %fprintf(fileID,'Xevent_sync_%s-%s EVENT_%s-%s_pos! EVENT_%s-%s_neg! EVENT!_%s-%s and2x4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});                                            
                % Este va, lo estoy sacando para probar el sincronizador
                %fprintf(fileID,'Xclk_sync_DFF_%s-%s sys_clk EVENT_S_%s-%s! clk_event_%s-%s clk_event_n_%s-%s sr_rw_sync vdd! dffsrx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});                  %clk=sys_clk %D=gnd! Q=EVENT! Qn=EVENTn! rs=clear_n ss=event_out!
                %prueba sincronizador
                fprintf(fileID,'Xclk_sync_DFF_%s-%s-1 sys_clk EVENT_S_%s-%s! clk_event_%s-%s-1 clk_event_n_%s-%s-1 sr_rw_sync vdd! dffsrx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});                  %clk=sys_clk %D=gnd! Q=EVENT! Qn=EVENTn! rs=clear_n ss=event_out!
                fprintf(fileID,'Xclk_sync_DFF_%s-%s-2 sys_clk clk_event_%s-%s-1 clk_event_%s-%s clk_event_n_%s-%s sr_rw_sync vdd! dffsrx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});                  %clk=sys_clk %D=gnd! Q=EVENT! Qn=EVENTn! rs=clear_n ss=event_out!

            else                                                         %LTspice piece of code
    %         	fprintf(fileID,'Aclk_sync EVENT_1! EVENT_2! clk_event_n VDD VDD VDD EVENT_S! 0 AND\n');    
    %           fprintf(fileID,'Aclk_sync_DFF EVENT_S! 0 sys_clk 0 resetDFFs clk_event_n clk_event 0 DFLOP Vhigh=1 Vlow=0 Trise=5n Tfall=5n Td=5n\n');                     %A2 D 0 clk pre clr Qn Q 0 DFLOP
    %         	fprintf(fileID,'Aevent_sync EVENT_1! EVENT_2! VDD VDD VDD VDD EVENT! 0 AND\n');    
            end
            
            
            fprintf(fileID,'\n');
            fprintf(fileID,'*Time-out counter\n');
            fprintf(fileID,'*************************************************************************\n');        
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                fprintf(fileID,'Xtime_out_%s-%s enable_read_pulse time_out_%s-%s clk_event_%s-%s ', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});
                for i=1:ceil(log2(max_pulses))
                    fprintf(fileID,'time_out_S_%s-%s%d-1 ', neural_layers{CPA_i}, neural_layers{CPA_i+1}, i);
                end                
                fprintf(fileID,'counter_%d_pulses\n', max_pulses);
                fprintf(fileID,'\n');  
                fprintf(fileID,'Xinv_TO_%s%s time_out_%s-%s time_out_n_%s-%s invx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});
            else
            end
            
            fprintf(fileID,'\n');
            fprintf(fileID,'*Read-write circuit synchronizer\n');
            fprintf(fileID,'*************************************************************************\n');        
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                fprintf(fileID,'X_WR_control_%s-%s write_ok!_%s-%s gnd! RW_%s-%s RWn_%s-%s resetDFFsn wr_user_n dffsrx4tf\n', neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1}, neural_layers{CPA_i}, neural_layers{CPA_i+1});                                                                          %D=write_ok! clk=sys_clk Q=EVENT! Qn=EVENTn! rs=clear_n ss=event_out!
            else
                %fprintf(fileID,'A_WR_control write_okn! 0 sys_clk wr_user resetDFFs RWn RW 0 DFLOP Vhigh=1 Vlow=0 Trise=5n Tfall=5n Td=5n\n');                                 %A2 D 0 clk pre clr Qn Q 0 DFLOP
                fprintf(fileID,'A_WR_control V_low2 0 write_ok! wr_user resetDFFs RWn RW 0 DFLOP Vhigh=1 Vlow=0 Trise=5n Tfall=5n Td=5n\n');                                    %A2 D 0 clk pre clr Qn Q 0 DFLOP
                fprintf(fileID,'Vlow2 V_low2 0 0\n');   
            end
            
        end
        
        fprintf(fileID,'Xlayer_sync ');
        for CPA_i=1:length(neural_layers)-1
            fprintf(fileID,'write_ok!_%s-%s ', neural_layers{CPA_i}, neural_layers{CPA_i+1});
        end
        fprintf(fileID,'write_ok! sys_clk vdd! resetDFFsn sync_W_OK_%d_layers\n\n', length(neural_layers)-1);

        
        %% Connections of the second layer
        fprintf(fileID,'\n');
            
        fprintf(fileID,'\n');
        fprintf(fileID,'*Read-write signals\n');
        fprintf(fileID,'*************************************************************************\n');        
        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
            fprintf(fileID,'Xinv_rw_reset resetDFFs resetDFFsn invx4tf\n');
            fprintf(fileID,'Xinv_rw_user wr_user wr_user_n invx4tf\n');    
        else
            %fprintf(fileID,'A_WR_control write_okn! 0 sys_clk wr_user resetDFFs RWn RW 0 DFLOP Vhigh=1 Vlow=0 Trise=5n Tfall=5n Td=5n\n');                                 %A2 D 0 clk pre clr Qn Q 0 DFLOP
            fprintf(fileID,'Vlow2 V_low2 0 0\n');
            %fprintf(fileID,'Ainv_write_ok write_ok! 0 0 0 0 write_okn! 0 0 BUF\n');
            fprintf(fileID,'A_WR_OK V_low2 0 sys_clk write_ok_pulse! resetDFFs write_okn! write_ok! 0 DFLOP Vhigh=1 Vlow=0 Trise=5n Tfall=5n Td=5n\n');                     %A2 D 0 clk pre clr Qn Q 0 DFLOP           
        end
       
        
        fprintf(fileID,'\n');            
        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
            fprintf(fileID,'*logic gates power suply\n');
            fprintf(fileID,'*************************************************************************\n');
            fprintf(fileID,'VVDD vdd! gnd! dc=1.2\n');
            fprintf(fileID,'VVSS vss! gnd! dc=0\n');
            fprintf(fileID,'\n');   
            fprintf(fileID,'*reser of DFFs\n');
            fprintf(fileID,'*************************************************************************\n');
            fprintf(fileID,'VresetDFFs resetDFFs gnd! dc=0 pulse ( 0 1.2 1u 10n 10n 1u %.12e )\n',9999);
            fprintf(fileID,'\n');    
            fprintf(fileID,'*System clock\n');
            fprintf(fileID,'*************************************************************************\n');
            if sim_WR==1
                fprintf(fileID,'Vclk sys_clk gnd! dc=0 pulse ( 0 1.2 0.5u 10n 10n %.5e %.5e )\n',sys_clk(1),sys_clk(2));
            else
                fprintf(fileID,'Vclk sys_clk gnd! dc=0 \n');
            end
            fprintf(fileID,'\n');                
        else
            fprintf(fileID,'*Comparator power suplies\n');
            fprintf(fileID,'*************************************************************************\n');
            fprintf(fileID,'VVDD VDD 0 5\n');
            fprintf(fileID,'VVSS VSS 0 0\n');
            fprintf(fileID,'\n');  
            fprintf(fileID,'*reser of DFFs\n');
            fprintf(fileID,'*************************************************************************\n');
            fprintf(fileID,'VresetDFFs resetDFFs 0 pulse ( 0 1.2 1u 10n 10n 1u 10 )\n');
            fprintf(fileID,'\n');    
            fprintf(fileID,'*System clock\n');
            fprintf(fileID,'*************************************************************************\n');
            fprintf(fileID,'Vclk sys_clk 0 pulse ( 0 1.2 0.5u 10n 10n 300u 600u )\n');
            fprintf(fileID,'\n');             
        end
        
        fprintf(fileID,'*Pulsed signal for memristor programming\n');
        fprintf(fileID,'*************************************************************************\n');       

        tON_write=1/wrFreq*DC-tr;
        tON_verify=1/wrFreq*(1-DC)*DC_v-tr;
        tOFF=1/wrFreq-(tON_write+tON_verify+2*tr);
        
        if sim_WR==1
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                  fprintf(fileID,'Vsetting-1 Vpulsedsignal Vpulsedsignal-1 dc=0 pulse ( 0 %f %.5g %.5g %.5g %.5g %.5g )\n',read_voltage, tdelay_gral, tr, tr, tON_verify, 1/wrFreq);    
                  fprintf(fileID,'Vsetting-2 Vpulsedsignal-1 gnd! dc=0 pulse ( 0 %f %.5g %.5g %.5g %.5g %.5g )\n',V_prog, tdelay_gral+tOFF/2+tON_verify+tr, tr, tr, tON_write, 1/wrFreq);   
                  fprintf(fileID,'Vsetting_en enable_read_pulse gnd! dc=0 pulse ( 0 %f %.5g %.5g %.5g %.5g %.5g )\n',VDD, tdelay_gral-tON_verify*0.05, tr, tr, tON_verify*1.1, 1/wrFreq);

%                 fprintf(fileID,'Vsetting-1 Vpulsedsignal Vpulsedsignal-1 dc=0 pulse ( 0 %f 0 1u 1u 10u 40u )\n',read_voltage);
%                 fprintf(fileID,'Vsetting-2 Vpulsedsignal-1 gnd! dc=0 pulse ( 0 %f 20u 1u 1u 10u 40u )\n',V_prog);        
%                 fprintf(fileID,'Vsetting_en enable_read_pulse gnd! dc=0 pulse ( 0 %f 0 0.1u 0.1u 11u 40u )\n',VDD);
                
            else
                fprintf(fileID,'Vsetting Vpulsedsignal 0 PULSE(0 %f 0 1u 1u 10u 20u)\n',V_prog);
            end
        else
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                fprintf(fileID,'Vsetting Vpulsedsignal gnd! dc=0\n');
                fprintf(fileID,'Vsetting_en enable_read_pulse gnd! dc=0\n');
            else
                fprintf(fileID,'Vsetting Vpulsedsignal 0 0\n');
            end
        end
        fprintf(fileID,'\n');
        
        fprintf(fileID,'*Read-Write signal\n');
        fprintf(fileID,'*************************************************************************\n'); 
        if sim_WR==1
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                fprintf(fileID,'V_WR WR_user gnd! dc=0 pulse (0 1 10u 10n 10n 1u %.12e)\n',9999);
            else
                fprintf(fileID,'V_WR WR_user 0 PULSE(0 1 10u 10n 10n 1u 10)\n');
            end
        else
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                fprintf(fileID,'V_WR WR_user gnd! dc=0\n');
            else
                fprintf(fileID,'V_WR WR_user 0 0\n');
            end
        end
        fprintf(fileID,'\n');
        fprintf(fileID,'*##########################################################################\n\n');

        fprintf(fileID,'*######################  .LIB / .INC statements  ##########################\n\n');
        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')

        else
            fprintf(fileID,'.lib UniversalOpamps2.sub\n');
            fprintf(fileID,'.model SW1 SW(Ron=1m Roff=100Meg Vt=.9 Vh=0\n');
            fprintf(fileID,'.model SW2 SW(Ron=1m Roff=100Meg Vt=-.1 Vh=0\n\n');            
        end
        fprintf(fileID,'*##########################################################################\n\n');

        fprintf(fileID,'*#########################  .PARAM statements  ############################\n\n');

        if size(Rs,2)==1
            fprintf(fileID,'.PARAM Rs=%f\n',Rs);    
        else
            fprintf(fileID,'.STEP PARAM Rs list ');    
            for i=1:size(Rs,2)
                fprintf(fileID,'%f ',Rs(i));
            end
            fprintf(fileID,'\n');
        end
        fprintf(fileID,'.PARAM Rs_PAD=%f\n',Rs_PAD);
        if calc_RE_latency==1
            if strcmpi(connections,'dualSide_connect') || strcmpi(connections,'dualSide_connect_outputs')
                n_extra=(neu_per_layer_sub(1)-neu_per_layer_sub(2))/2;
            else
                n_extra=(neu_per_layer_sub(1)-neu_per_layer_sub(2));
            end
        else
            n_extra=1;            
        end
        fprintf(fileID,'.PARAM n_extra=%d\n',n_extra);
        fprintf(fileID,'.PARAM Rcs=%f\n',R_cs);
        fprintf(fileID,'.PARAM R_shunt=%f\n',R_shunt);
        fprintf(fileID,'.PARAM HZ=%s\n',R_off_value);
        fprintf(fileID,'.PARAM R_sense=%f\n',R_sense);
        fprintf(fileID,'.PARAM R_in_TIA=%f\n',R_in_TIA);
        if strcmpi(use_caps,'yes')
            fprintf(fileID,'.PARAM C_interline=%e\n',C_interline);  
            fprintf(fileID,'.PARAM C_line2gnd=%e\n',C_line2gnd);  
            fprintf(fileID,'.PARAM C_memdiode=%e\n',C_memdiode);   
        end
        fprintf(fileID,'*##########################################################################\n\n');

        fprintf(fileID,'*#########################  .SIM statements  ############################\n\n');
        if strcmp(simEngine,'HSPICE')
            fprintf(fileID,'.tran %.12e %.12e 0\n\n',max_time_step, simulation_time);
        elseif strcmp(simEngine,'FineSim')
            if strcmpi(noise_transient,'yes')
                fprintf(fileID,'.tran %.12e %.12e 0 noisefmax=%.2e noisefmin=%.2e noisescale=%.2e\n\n',max_time_step, simulation_time, noisefmax, noisefmin, noisescale);
            else
                fprintf(fileID,'.tran %.12e %.12e 0\n\n',max_time_step, simulation_time);
            end
        else
            fprintf(fileID,'.tran 0 %.12e 0 %.12e\n\n',simulation_time,max_time_step);
        end
        fprintf(fileID,'*##########################################################################\n\n');

        % .PRINT statements
        fprintf(fileID,'*#########################  .PRINT statements  ############################\n\n');
        
        if strcmp(analysis,'transient')
            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                if strcmp(output_fmt,'binary')
                    save_statement='.PROBE';
                else
                    save_statement='.PRINT';
                end
                %fprintf(fileID,'.PRINT TRAN ');
            else
                save_statement='.PRINT';
                fprintf(fileID,'.PRINT TRAN time ');
            end
            
            if strcmpi(debug_mode,'yes')
                fprintf(fileID,'%s TRAN V(*)\n',save_statement);
                fprintf(fileID,'%s TRAN I(*)\n',save_statement);                
            else
                for CPA_i=1:length(neural_layers)-1
                    for pos_neg_i=1:2
                        for partitions_i=1:partitions_N(CPA_i)
                            for partitions_ii=1:partitions_N(CPA_i+1)
                                for i=1:neu_per_layer_sub(CPA_i)
                                    for ii=1:neu_per_layer_sub(CPA_i+1)
                                        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                            if save_all==1 || strcmpi(save_h,'yes')
                                                fprintf(fileID,'%s TRAN V(xmemd_network_%s-%s_P%d-%d_%s.xmemdr%d-%dc%d-%d.h)\n',save_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii,ii,i);
                                            end
                                            %I(mmed_network:Rcsa1-1b1-1_1)
                                            if save_all==1  || strcmpi(save_imemd,'yes')
                                                fprintf(fileID,'%s TRAN I(xmemd_network_%s-%s_P%d-%d_%s.Rcs_r%d-%dc%d-%d_1)\n',save_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii,ii,i);
                                            end
                                            %V(mmed_network:a1-1_p)
                                            if save_all==1  || strcmpi(save_vpn,'yes')
                                                fprintf(fileID,'%s TRAN V(xmemd_network_%s-%s_P%d-%d_%s.r%d-%d_p)\n',save_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii);
                                                %V(mmed_network:b1-1_n)
                                                fprintf(fileID,'%s TRAN V(xmemd_network_%s-%s_P%d-%d_%s.c%d-%d_n)\n',save_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},ii,i);
                                            end
                                            %fprintf(fileID,'\n');
                                        else
                                            fprintf(fileID,'V(memd_network_P%d_%s:memdr%d-%dc%d-%d:h) ',partitions_i, pos_neg{pos_neg_i},i,ii,ii,i);
                                            fprintf(fileID,'I(memd_network_P%d_%s:Rcsa%d-%db%d-%d_1) ',partitions_i, pos_neg{pos_neg_i},i,ii,ii,i);
                                            fprintf(fileID,'V(memd_network_P%d_%s:a%d-%d_p) ',partitions_i, pos_neg{pos_neg_i},i,ii);                    
                                            fprintf(fileID,'V(memd_network_P%d_%s:b%d-%d_n) ',partitions_i, pos_neg{pos_neg_i},ii,i);
                                        end
                                    end
                                end
                                
                                %I(Rneuron<layer>1)
                                if strcmpi(save_inner_neuron,'yes') || save_all==1
                                    for i=1:neu_per_layer_sub(CPA_i+1) 
                                        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                            %fprintf(fileID,'%s TRAN I(Rneuron%s%d_P%d-%d_%s) \n',save_statement,neural_layers{CPA_i+1},i,partitions_i, partitions_ii, pos_neg{pos_neg_i});
                                            fprintf(fileID,'%s TRAN i__Rneuron%s%d_P%d-%d_%s__=PAR(''1*ISUB(Xselector_RW_%s%s%d-1-P%d-%d_%s.Neuron)'') \n',save_statement,neural_layers{CPA_i+1},i,partitions_i, partitions_ii, pos_neg{pos_neg_i},neural_layers{CPA_i+1},neural_layers{CPA_i},i,partitions_i, partitions_ii, pos_neg{pos_neg_i});
                                            %fprintf(fileID,'%s TRAN ISUB(Xselector_RW_%s%s%d-1-P%d-%d_%s.Neuron) \n',save_statement,neural_layers{CPA_i+1},neural_layers{CPA_i},i,partitions_i, partitions_ii, pos_neg{pos_neg_i});
                                        else
                                            fprintf(fileID,'I(Rneuron%s%d_P%d_%s) ',neural_layers{2},i,partitions_i, pos_neg{pos_neg_i});
                                        end
                                        fprintf(fileID,'\n');
                                    end
                                end
                                
                                %V(Rneuron<layer>1)
                                if pos_neg_i==2
                                    if  CPA_i==length(neural_layers)-1
                                        if partitions_i==partitions_N(CPA_i)
                                            for i=1:neu_per_layer_sub(CPA_i+1) 
                                                if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                                    fprintf(fileID,'%s TRAN V(Neuron%s%d_P%d-1) \n',save_statement,neural_layers{CPA_i+1},i, partitions_ii);
                                                else
                                                end
                                                fprintf(fileID,'\n');
                                            end
                                        end
                                    end   
                                end
                                
                                
    
                                if save_all==1 || strcmpi(save_write_signals,'yes')
                                    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                        fprintf(fileID,'%s TRAN V(enable_write_%s-%s_P%d-%d_%s)\n',save_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}); 
                                        fprintf(fileID,'%s TRAN V(comparator_%s-%s_P%d-%d_%s)\n',save_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                                        fprintf(fileID,'%s TRAN V(comparator_%s-%s_P%d-%d_%s-1)\n',save_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii , pos_neg{pos_neg_i});
                                        fprintf(fileID,'%s TRAN V(EVENT_%s-%s_P%d-%d_%s!)\n',save_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i});
                                    else
                                    end
                                end
                            end
                        end
                        if save_all==1 || strcmpi(save_write_signals,'yes')
                            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                fprintf(fileID,'%s TRAN V(EVENT_%s-%s_%s!)\n',save_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i});
                            else
                            end
                        end
                    end
                end
                if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                    fprintf(fileID,'\n%s TRAN V(write_ok!)\n',save_statement);
                    fprintf(fileID,'%s TRAN V(sys_clk)\n',save_statement);
                    fprintf(fileID,'%s TRAN V(RW)\n',save_statement);   
                    fprintf(fileID,'%s TRAN V(eventsA)\n',save_statement);
                    fprintf(fileID,'%s TRAN V(Vpulsedsignal)\n',save_statement);
                    fprintf(fileID,'%s TRAN V(enable_read_pulse)\n',save_statement);

                else
                    %fprintf(fileID,'\n.PRINT TRAN V(write_ok!)\n');
                    %fprintf(fileID,'.PRINT TRAN V(sys_clk)\n');
                    %fprintf(fileID,'.PRINT TRAN V(enable_write_1)\n');  
                    %fprintf(fileID,'.PRINT TRAN V(enable_write_2)\n');  
                    %fprintf(fileID,'.PRINT TRAN V(RW)\n');   
                    %fprintf(fileID,'.PRINT TRAN V(comparator_1)\n');      
                    %fprintf(fileID,'.PRINT TRAN V(comparator_2)\n');                  
                    %fprintf(fileID,'.PRINT TRAN V(eventsA)\n');      
                    %fprintf(fileID,'.PRINT TRAN V(event_1!)\n');
                    %fprintf(fileID,'.PRINT TRAN V(event_2!)\n');                
                end
            end

        elseif strcmp(analysis,'o_point')
            for CPA_i=1:length(neural_layers)-1
                for pos_neg_i=1:2
                    for partitions_i=1:partitions_N(CPA_i)
                        for partitions_ii=1:partitions_N(CPA_i+1)
                            %I(mmed_network:Rcsa1-1b1-1_1)
                            for i=1:neu_per_layer_sub(CPA_i)
                                for ii=1:neu_per_layer_sub(CPA_i+1)
                                    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                        fprintf(fileID,'%s TRAN I(xmemd_network_%s-%s_P%d-%d_%s.Rcsr%d-%dc%d-%d_1) FIND I(xmemd_network_%s-%s_P%d-%d_%s.Rcsr%d-%dc%d-%d_1) AT=%.12e\n',meas_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii,ii,i,neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii,ii,i,simulation_time);                
                                    else
                                        fprintf(fileID,'I(mmed_network_P%d_%s:Rcsa%d-%db%d-%d_1) ', partitions_i, pos_neg{pos_neg_i},i,ii,ii,i);
                                    end
                                end
                            end

                            %V(mmed_network:b1-1_pn)
                            for i=1:neu_per_layer_sub(CPA_i)
                                for ii=1:neu_per_layer_sub(CPA_i+1)
                                    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                        fprintf(fileID,'%s TRAN V(xmemd_network_%s-%s_P%d-%d_%s.r%d-c%d_pn) FIND PAR(''V(xmemd_network_%s-%s_P%d-%d_%s.r%d-%d_p)-V(xmemd_network_%s-%s_P%d-%d_%s.c%d-%d_n)'') AT=%.12e\n',meas_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii,neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},ii,i,simulation_time);
                                    else
                                        fprintf(fileID,'V(mmed_network_P%d_%s:b%d-%d_n) ', partitions_i, pos_neg{pos_neg_i},i,ii);
                                    end
                                end
                            end   
                        end
                    end
                end
            end
        end
        if (strcmp(analysis,'transient') && sim_WR==1) || (strcmp(analysis,'transient') && strcmpi(save_h,'yes'))
            for CPA_i=1:length(neural_layers)-1
                for pos_neg_i=1:2
                    for partitions_i=1:partitions_N(CPA_i)
                        for partitions_ii=1:partitions_N(CPA_i+1)
                            %V(mmed_network:b1-1_pn)
                            for i=1:neu_per_layer_sub(CPA_i)
                                for ii=1:neu_per_layer_sub(CPA_i+1)
                                    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')  
                                        fprintf(fileID,'%s TRAN V(xmemd_network_%s-%s_P%d-%d_%s.xmemdr%d-%dc%d-%d.h) FIND V(xmemd_network_%s-%s_P%d-%d_%s.xmemdr%d-%dc%d-%d.h) AT=%.12e\n',meas_statement, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii,ii,i, neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii,ii,i,simulation_time);
                                    else
                                        fprintf(fileID,'V(mmed_network_P%d_%s:b%d-%d_n) ', partitions_i, pos_neg{pos_neg_i},i,ii);
                                    end
                                end
                            end 
                        end
                    end
                end      
            end
        end
        if strcmp(analysis,'transient') && calc_WR_margin==1
            for pos_neg_i=1:2
                for partitions_i=1:partitions_N
                    %V(mmed_network:b1-1_pn)
                    for i=1:neu_per_layer_sub(1)
                        for ii=1:neu_per_layer_sub(2)
                            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')  
                                fprintf(fileID,'%s TRAN V(xmemd_network_P%d_%s.a%d-b%d_pn) MAX PAR(''V(xmemd_network_P%d_%s.a%d-%d_p)-V(xmemd_network_P%d_%s.b%d-%d_n)'')\n',meas_statement,partitions_i, pos_neg{pos_neg_i},i,ii,partitions_i, pos_neg{pos_neg_i},i,ii,partitions_i, pos_neg{pos_neg_i},ii,i);
                            else
                                fprintf(fileID,'V(mmed_network_P%d_%s:b%d-%d_n) ', partitions_i, pos_neg{pos_neg_i},i,ii);
                            end
                        end
                    end                 
                end
            end        
        end
        if strcmp(analysis,'transient') && strcmpi(calc_pwr,'yes')
            for power_meas_statement_idx=1:length(power_meas_statement)
                if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')  
                    fprintf(fileID,'%s\n',power_meas_statement{power_meas_statement_idx,1});
                else
                    fprintf(fileID,'V(mmed_network_P%d_%s:b%d-%d_n) ', partitions_i, pos_neg{pos_neg_i},i,ii);
                end
            end        
        end
        if strcmp(analysis,'transient') && calc_RE_latency==1
%             for pos_neg_i=1:2
%                 for partitions_i=1:partitions_N
%                     %I(Rneuronb1)
%                     for ii=1:neu_per_layer_sub(2)
%                         if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
%                             fprintf(fileID,'%s TRAN I_Rneuron%s%d_P%d_%s_final FIND I(Rneuron%s%d_P%d_%s) AT=%.12e\n',meas_statement,neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},simulation_time);
%                             fprintf(fileID,'%s TRAN I_Rneuron%s%d_P%d_%s_settling_max WHEN I(Rneuron%s%d_P%d_%s)=''I_Rneuron%s%d_P%d_%s_final*1.05'' CROSS=LAST\n',meas_statement,neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i});
%                             fprintf(fileID,'%s TRAN I_Rneuron%s%d_P%d_%s_settling_min WHEN I(Rneuron%s%d_P%d_%s)=''I_Rneuron%s%d_P%d_%s_final*0.95'' CROSS=LAST\n',meas_statement,neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i});
%                         else
%                             fprintf(fileID,'I(Rneuron%s%d_P%d_%s) ',neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i});
%                         end
%                     end
%                 end
%             end
            
            for pos_neg_i=1:2
                for partitions_i=1:partitions_N
                %I(mmed_network:Rcsa1-1b1-1_1)
                    for i=1:neu_per_layer_sub(1)
                        for ii=1:neu_per_layer_sub(2)
                            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                if i==neu_per_layer_sub(1)
                                    fprintf(fileID,'%s TRAN I_xmemd_network_P%d_%s.Rcsa%d_%db%d_%d_1_final FIND I(xmemd_network_P%d_%s.Rcsa%d-%db%d-%d_1) AT=%.12e\n',meas_statement, partitions_i, pos_neg{pos_neg_i},i,ii,ii,i,partitions_i, pos_neg{pos_neg_i},i,ii,ii,i,simulation_time);                
                                    fprintf(fileID,'%s TRAN I_xmemd_network_P%d_%s.Rcsa%d_%db%d_%d_1_settling_max WHEN I(xmemd_network_P%d_%s.Rcsa%d-%db%d-%d_1)=''I_xmemd_network_P%d_%s.Rcsa%d_%db%d_%d_1_final*%f'' CROSS=LAST\n',meas_statement, partitions_i, pos_neg{pos_neg_i},i,ii,ii,i,partitions_i, pos_neg{pos_neg_i},i,ii,ii,i, partitions_i, pos_neg{pos_neg_i},i,ii,ii,i,1+error_band);
                                    fprintf(fileID,'%s TRAN I_xmemd_network_P%d_%s.Rcsa%d_%db%d_%d_1_settling_min WHEN I(xmemd_network_P%d_%s.Rcsa%d-%db%d-%d_1)=''I_xmemd_network_P%d_%s.Rcsa%d_%db%d_%d_1_final*%f'' CROSS=LAST\n',meas_statement, partitions_i, pos_neg{pos_neg_i},i,ii,ii,i,partitions_i, pos_neg{pos_neg_i},i,ii,ii,i, partitions_i, pos_neg{pos_neg_i},i,ii,ii,i,1-error_band);
                                end
                                
                                if i==1
                                    fprintf(fileID,'%s TRAN I_Rneuron%s%d_P%d_%s_final FIND I(Rneuron%s%d_P%d_%s) AT=%.12e\n',meas_statement,neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},simulation_time);
                                    fprintf(fileID,'%s TRAN I_Rneuron%s%d_P%d_%s_settling_max WHEN I(Rneuron%s%d_P%d_%s)=''I_Rneuron%s%d_P%d_%s_final*%f'' CROSS=LAST\n',meas_statement,neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},1+error_band);
                                    fprintf(fileID,'%s TRAN I_Rneuron%s%d_P%d_%s_settling_min WHEN I(Rneuron%s%d_P%d_%s)=''I_Rneuron%s%d_P%d_%s_final*%f'' CROSS=LAST\n',meas_statement,neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i},1-error_band);
                                end
                            
                            
                            else
                                fprintf(fileID,'I(mmed_network_P%d_%s:Rcsa%d-%db%d-%d_1) ', partitions_i, pos_neg{pos_neg_i},i,ii,ii,i);
                            end
                        end
                    end                 
                end
            end
            
        end
        
        fprintf(fileID,'\n\n');
        fprintf(fileID,'*########################################################################\n\n');

        fprintf(fileID,'.end');

        fclose(fileID);
        tElapsed_netlist=toc(tStart_netlist);
        
        fprintf('done!\n');
        if guiStatus_data
            str_to_display{end}=sprintf("%sCreating circuit netlist for a Neural Network of %d by %d...done! (%s)",str_header, neu_per_layer(1), neu_per_layer(2),duration([0, 0, tElapsed_netlist]));
            figure(status_gui);
            status_guiObject.log.String=str_to_display;  
            status_guiObject.log.Value=length(str_to_display)+1;
            drawnow
        end
        
    %% Netlist simulation
    % The following piece of code performs the circuit simulation depenging
    % on the chosen simEngine. It also process and stores the resuls
    
        if strcmp(simulate_netlist,'yes')
            tStart_sim=tic;            
            sim_complete=0;
            while ~sim_complete
                if guiStatus_data
                    str_to_display{end+1}=sprintf("%sSimulating circuit netlist (%.3e sec.)...",str_header,simulation_time);
                    figure(status_gui);
                    status_guiObject.log.String=str_to_display;
                    status_guiObject.log.Value=length(str_to_display)+1;
                    drawnow
                end
                
                tStart_curr_sim=tic;

                if exist(results_folder,'dir')==0
                    mkdir(results_folder);
                end                
                if ispc
                    exe_file="C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe";
                    system_command = strcat('"',exe_file,'" -run',{' '},runMode,{' '},fullfile(sim_folder,netlist_name));
                    [status,cmdout]=system(system_command,'-echo');
                elseif isunix
                    if strcmp(simEngine,'HSPICE')
                        %unix('module load HSPICE');
                        exe_file=strcat(hspice_root_dir,'bin',filesep,'hspice');
                        system_command = strcat('"',exe_file,'" -i',{' '},fullfile(sim_folder,netlist_name),{' '},'-n 1',{' '},'-o',{' '},fullfile(sim_folder,strrep(netlist_name,'.net','')));
                        system_command=system_command{1};
                        [status,cmdout]=unix(system_command,'-echo');    
                    elseif strcmp(simEngine,'FineSim')
                        exe_file=strcat(finesim_dir,'bin',filesep,'finesim');
                        system_command = strcat('"',exe_file,'" -i',{' '},fullfile(sim_folder,netlist_name),{' '},'-np',{' '},num2str(nProc),{' '},'-o',{' '},fullfile(sim_folder,strrep(netlist_name,'.net','')));
                        system_command=system_command{1};
                        [status,cmdout]=unix(system_command,'-echo');                         
                    else
                        exe_file="C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe";
                        system_command = strcat('"',exe_file,'" -run',{' '},runMode,{' '},fullfile(sim_folder,netlist_name));
                        [status,cmdout]=system(system_command,'-echo');                        
                    end
                end

                if sim_WR==0
                    tElapsed_sim_read=toc(tStart_sim);
                else
                    tElapsed_curr_sim=toc(tStart_curr_sim);
                end
                
                if strcmp(simEngine,'FineSim') || strcmp(simEngine,'HSPICE')
                % This part of the codes process the simulated data in case of having chosen the  FineSim simulation engine
                    if strcmp(simEngine,'FineSim')
                        if ~strcmp(sim_folder,results_folder)
                            if exist(fullfile(results_folder,'figures'),'dir')
                                rmdir(fullfile(results_folder,'figures'),'s');   
                            end   
                            if exist(fullfile(results_folder,'input_sources'),'dir')
                                rmdir(fullfile(results_folder,'input_sources'),'s'); 
                            end

                            if exist(fullfile(results_folder,'signals'),'dir')
                                rmdir(fullfile(results_folder,'signals'),'s'); 
                            end    

                            if exist(fullfile(results_folder,'defaults'),'dir')
                                rmdir(fullfile(results_folder,'defaults'),'s'); 
                            end   

                            delete(fullfile(results_folder,'*'));

                            movefile(fullfile(sim_folder,strrep(netlist_name,'.net','*')),results_folder);
                            movefile(fullfile(sim_folder,'input_sources'),results_folder);
                            copyfile(fullfile('..','config_files'),results_folder);
                        end   
                    elseif strcmp(simEngine,'HSPICE')
                        if ~strcmp(sim_folder,results_folder)
                            if exist(fullfile(results_folder,strrep(netlist_name,'.net','.pvadir')),'dir')
                                rmdir(fullfile(results_folder,strcat(strrep(netlist_name,'.net','*'),'.pvadir')),'s');   
                            end
                            if exist(fullfile(results_folder,'figures'),'dir')
                                rmdir(fullfile(results_folder,'figures'),'s');   
                            end 
                            if exist(fullfile(results_folder,'input_sources'),'dir')
                                rmdir(fullfile(results_folder,'input_sources'),'s'); 
                            end

                            if exist(fullfile(results_folder,'signals'),'dir')
                                rmdir(fullfile(results_folder,'signals'),'s'); 
                            end   
                            if exist(fullfile(results_folder,'defaults'),'dir')
                                rmdir(fullfile(results_folder,'defaults'),'s'); 
                            end   

                            delete(fullfile(results_folder,'*'));
                            if exist(fullfile(sim_folder,strrep(netlist_name,'.net','.printtr0')),'file')
                                delete(fullfile(sim_folder,strrep(netlist_name,'.net','.printtr0')));
                            end

                            movefile(fullfile(sim_folder,strrep(netlist_name,'.net','*')),results_folder);
                            movefile(fullfile(sim_folder,'input_sources'),results_folder);
                            copyfile(fullfile('..','config_files'),results_folder);

                            if exist(fullfile(sim_folder,strrep(netlist_name,'.net','.pvadir')),'dir')
                                rmdir(fullfile(sim_folder,strcat(strrep(netlist_name,'.net','*'),'.pvadir')),'s');   
                            end
                        end
                    end
                    
                    if strcmp(analysis,'o_point')
                            raw_data_filename=strcat('NW_',num2str(neu_per_layer(1)),'by',num2str(neu_per_layer(2)),'_closed_loop.mt0');
                            map=read_opData_part(fullfile(results_folder,raw_data_filename),partitions_N, results_folder,'networkSize',neu_per_layer);
                            map.folder=results_folder;
                            sim_complete=1;
                    else                         
                        if strcmpi(output_fmt,'ascii')
                            if strcmp(simEngine,'FineSim')
                                raw_data_filename=strcat('NW_',num2str(neu_per_layer(1)),'by',num2str(neu_per_layer(2)),'_closed_loop.pt0');                    
                                RAW_DATA=read_FineSimresults(fullfile(results_folder,raw_data_filename));
                                RAW_DATA.neu_per_layer=neu_per_layer;  
                                RAW_DATA.simEngine='FineSim';
                                save(fullfile(results_folder,'RAW_DATA.mat'),'RAW_DATA','W_matrix');  
                                sim_complete=1;                               
                            elseif strcmp(simEngine,'HSPICE')
                                raw_data_filename=strcat('NW_',num2str(neu_per_layer(1)),'by',num2str(neu_per_layer(2)),'_closed_loop.tr1');
                                RAW_DATA=read_HSPICEresults(fullfile(results_folder,raw_data_filename));
                                RAW_DATA.neu_per_layer=neu_per_layer;  
                                RAW_DATA.simEngine='HSPICE';
                                save(fullfile(results_folder,'RAW_DATA.mat'),'RAW_DATA','W_matrix');        
                            end

                        elseif strcmpi(output_fmt,'binary')

                            read_FSDB(results_folder, simEngine, wv_dir, netlist_name, neural_layers, sim_WR, remove_RAW_after_sim);
                                
                            if sim_WR==1

                                run(fullfile(results_folder,strrep(netlist_name,'.net','.m')))
                                
                                WR_OK_pos=find(v__write_ok___>0.6);
                                if isempty(WR_OK_pos)
                                    repeat_write_sim=1;
                                else
                                    if TIME(1,max(WR_OK_pos))<1e-6
                                        repeat_write_sim=1;
                                    else
                                        repeat_write_sim=0;
                                    end
                                end
                                
                                if repeat_write_sim==0
                                    if calc_WR_margin==0
                                        
                                        raw_data_filename=strcat(dir_n_files.NW_files_id,'.mt0');
                                        map=read_hData_MLP(fullfile(results_folder,raw_data_filename),partitions_N, neural_layers, results_folder,'h_DATA','networkSize',neu_per_layer);
                                        map.folder=results_folder;

                                        WR_OK_time=TIME(1,max(WR_OK_pos));
                                        map.write_time=WR_OK_time;
                                        map.simulation_time=simulation_time;
                                        map.W_matrix_init_pos=W_matrix_init_pos;
                                        map.W_matrix_init_neg=W_matrix_init_neg;

                                        str_to_display{end}=sprintf("%sSimulating circuit netlist (%.3e sec.)...done! (sim time: %.3g sec, run time: %s)",str_header,simulation_time,WR_OK_time,duration([0, 0, tElapsed_curr_sim]));
                                        figure(status_gui);
                                        status_guiObject.log.String=str_to_display;
                                        status_guiObject.log.Value=length(str_to_display)+1;
                                        drawnow

                                        save(fullfile(results_folder,'h_matrix.mat'),'map');
                                        sim_complete=1;
                                        
                                    elseif calc_WR_margin==1
                                        
                                        raw_data_filename=strcat('NW_',num2str(neu_per_layer(1)),'by',num2str(neu_per_layer(2)),'_closed_loop.mt0');
                                        map=read_hData_part(fullfile(results_folder,raw_data_filename),partitions_N, results_folder,'vPN_DATA','networkSize',neu_per_layer);
                                        map.folder=results_folder;

                                        WR_OK_time=TIME(1,max(WR_OK_pos));
                                        map.write_time=WR_OK_time;
                                        map.simulation_time=simulation_time;

                                        str_to_display{end}=sprintf("%sSimulating circuit netlist (%.3e sec.)...done! (sim time: %.3g sec, run time: %s)",str_header,simulation_time,WR_OK_time,duration([0, 0, tElapsed_curr_sim]));
                                        figure(status_gui);
                                        status_guiObject.log.String=str_to_display;
                                        status_guiObject.log.Value=length(str_to_display)+1;
                                        drawnow

                                        save(fullfile(results_folder,'WR_margins.mat'),'map');
                                        sim_complete=1;
                                    end
                                    
                                else
                                    sim_complete=0;
                                    str_to_display{end}=sprintf("%sSimulating circuit netlist (%.3e sec.)...failed. Increasing sim. time (run time: %s)",str_header,simulation_time,duration([0, 0, tElapsed_curr_sim]));
                                    figure(status_gui);
                                    status_guiObject.log.String=str_to_display;
                                    status_guiObject.log.Value=length(str_to_display)+1;
                                    drawnow

                                    simulation_time_aux=simulation_time*1.5;
                                    SearchString=sprintf('.tran %.12e %.12e 0',max_time_step,simulation_time);
                                    ReplaceString=sprintf('.tran %.12e %.12e 0',max_time_step,simulation_time_aux);
                                    func_replace_string(fullfile(results_folder,netlist_name), fullfile(sim_folder,netlist_name), SearchString, ReplaceString)
                                    simulation_time=simulation_time_aux;  
                                end
                            else
                                sim_complete=1;
                                str_to_display{end}=sprintf("%sSimulating circuit netlist (%.3e sec.)...done! (run time: %s)",str_header,simulation_time,duration([0, 0, tElapsed_sim_read]));
                                figure(status_gui);
                                status_guiObject.log.String=str_to_display;
                                status_guiObject.log.Value=length(str_to_display)+1;
                                drawnow
                                inference_vars.sim_time=tElapsed_sim_read;
                                inference_vars.SPICE_sim_tim=simulation_time;
                                inference_vars.SPICE_RUNLVL=RUNLVL;
                                inference_vars.W_pos_init=W_matrix_init_pos;
                                inference_vars.W_neg_init=W_matrix_init_neg;
 
                                inference_vars.W_pos_init_fresh=CPA_settings.W_init_pos_fresh;
                                inference_vars.W_neg_init_fresh=CPA_settings.W_init_neg_fresh;     
                                
                                save(fullfile(results_folder,'inference_status.mat'),'inference_vars');
                                
                                if (strcmp(analysis,'transient') && strcmpi(save_h,'yes'))
                                    raw_data_filename=strcat(dir_n_files.NW_files_id,'.mt0');
                                    map=read_hData_MLP(fullfile(results_folder,raw_data_filename),partitions_N, neural_layers, results_folder,'h_DATA','networkSize',neu_per_layer);
                                    map.folder=results_folder;

%                                     WR_OK_time=TIME(1,max(WR_OK_pos));
%                                     map.write_time=WR_OK_time;
                                    map.simulation_time=simulation_time;
                                    map.W_matrix_init_pos=W_matrix_init_pos;
                                    map.W_matrix_init_neg=W_matrix_init_neg;

                                    save(fullfile(results_folder,'h_matrix.mat'),'map');
                                end
                                
                                if strcmpi(calc_pwr,'yes')
                                    if exist(fullfile(results_folder,strrep(netlist_name,'.net','.mt0')),'file')
                                        pwr_data=importdata(fullfile(results_folder,strrep(netlist_name,'.net','.mt0')));
                                        RU_idx=[];
                                        RC_idx=[];
                                        RL_idx=[];
                                        MD_idx=[];

                                        for col_idx=1:length(pwr_data.colheaders)
                                            if strfind(pwr_data.colheaders{1,col_idx},'.ru')
                                                RU_idx=horzcat(RU_idx,col_idx);
                                            end
                                            if strfind(pwr_data.colheaders{1,col_idx},'.rc')
                                                RC_idx=horzcat(RC_idx,col_idx);
                                            end
                                            if strfind(pwr_data.colheaders{1,col_idx},'.rl')
                                                RL_idx=horzcat(RL_idx,col_idx);
                                            end 
                                            if strfind(pwr_data.colheaders{1,col_idx},'.md')
                                                MD_idx=horzcat(MD_idx,col_idx);
                                            end                                                 
                                        end

                                        MD_pwr=sum(pwr_data.data(1,MD_idx));
                                        R_pwr=sum(pwr_data.data(1,RU_idx))+sum(pwr_data.data(1,RC_idx))+sum(pwr_data.data(1,RL_idx));

                                        % esto es para calcular la cantidad
                                        % de operaciones por joule...en
                                        % desarrollo
                                        total_num_dispo=0;
                                        for cnt_num_dispos_i=1:length(neu_per_layer)-1
                                            total_num_dispo=total_num_dispo+neu_per_layer(cnt_num_dispos_i)*neu_per_layer(cnt_num_dispos_i+1)*2;
                                        end
                                        op_per_joule=((simulation_time*input_vector_freq)*total_num_dispo)/(simulation_time*(MD_pwr+R_pwr));
                                        %
                                        
                                        save(fullfile(results_folder,'calc_pwr.mat'),'MD_pwr','R_pwr','pwr_data','op_per_joule');
                                    end
                                elseif calc_RE_latency==1
                                    if exist(fullfile(results_folder,strrep(netlist_name,'.net','.mt0')),'file')
                                        latency_data=importdata(fullfile(results_folder,strrep(netlist_name,'.net','.mt0')));
                                        
                                        RE_latency_file=fopen(fullfile(results_folder,strrep(netlist_name,'.net','.mt0')));
                                        RE_latency_line=fgetl(RE_latency_file);
                                        while ~feof(RE_latency_file)
                                            if ~isempty(strfind(RE_latency_line,'settling_max')) && ~isempty(strfind(RE_latency_line,'settling_min')) && ~isempty(strfind(RE_latency_line,'final'))
                                                headers=split(RE_latency_line);
                                                RE_latency_line=fgetl(RE_latency_file);
                                                data=split(RE_latency_line);
                                            else
                                                RE_latency_line=fgetl(RE_latency_file);
                                            end
                                        end
                                        for header_idx=1:size(headers,1)
                                            
                                            for pos_neg_i=1:2
                                                for partitions_i=1:partitions_N
                                                %I(mmed_network:Rcsa1-1b1-1_1)
                                                    for i=1:neu_per_layer_sub(1)
                                                        for ii=1:neu_per_layer_sub(2)
                                                            if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                                                if i==neu_per_layer_sub(1)
                                                                    string2search_final=sprintf('I_xmemd_network_P%d_%s.Rcsa%d_%db%d_%d_1_final', partitions_i, pos_neg{pos_neg_i},i,ii,ii,i);                
                                                                    string2search_settle_min=sprintf('I_xmemd_network_P%d_%s.Rcsa%d_%db%d_%d_1_settling_max', partitions_i, pos_neg{pos_neg_i},i,ii,ii,i);
                                                                    string2search_settle_max=sprintf('I_xmemd_network_P%d_%s.Rcsa%d_%db%d_%d_1_settling_min', partitions_i, pos_neg{pos_neg_i},i,ii,ii,i);
                                                                    str2search_neuron_final=sprintf('I_Rneuron%s%d_P%d_%s_final',neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i});
                                                                    str2search_neuron_settle_max=sprintf('I_Rneuron%s%d_P%d_%s_settling_max',neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i});
                                                                    str2search_neuron_settle_min=sprintf('I_Rneuron%s%d_P%d_%s_settling_min',neural_layers{2},ii,partitions_i, pos_neg{pos_neg_i});

                                                                    if strcmpi(headers{header_idx,1},string2search_final)
                                                                    elseif strcmpi(headers{header_idx,1},string2search_settle_min)
                                                                        if ~strcmpi(data{header_idx,1},'error')
                                                                            settle_min(i,ii,partitions_i,pos_neg_i)=str2num(data{header_idx,1});
                                                                        else
                                                                            settle_min(i,ii,partitions_i,pos_neg_i)=0;
                                                                        end
                                                                    elseif strcmpi(headers{header_idx,1},str2search_neuron_settle_min)
                                                                        if ~strcmpi(data{header_idx,1},'error')
                                                                            settle_neuron_min(ii,partitions_i,pos_neg_i)=str2num(data{header_idx,1});
                                                                        else
                                                                            settle_neuron_min(ii,partitions_i,pos_neg_i)=0;
                                                                        end
                                                                    elseif strcmpi(headers{header_idx,1},string2search_settle_max)
                                                                        if ~strcmpi(data{header_idx,1},'error')
                                                                            settle_max(i,ii,partitions_i,pos_neg_i)=str2num(data{header_idx,1});
                                                                        else
                                                                            settle_max(i,ii,partitions_i,pos_neg_i)=0;
                                                                        end
                                                                    elseif strcmpi(headers{header_idx,1},str2search_neuron_settle_max)
                                                                        if ~strcmpi(data{header_idx,1},'error')
                                                                            settle_neuron_max(ii,partitions_i,pos_neg_i)=str2num(data{header_idx,1});
                                                                        else
                                                                            settle_neuron_max(ii,partitions_i,pos_neg_i)=0;
                                                                        end
                                                                    end
                                                                end
                                                            else
                                                                fprintf(fileID,'I(mmed_network_P%d_%s:Rcsa%d-%db%d-%d_1) ', partitions_i, pos_neg{pos_neg_i},i,ii,ii,i);
                                                            end
                                                        end
                                                    end                 
                                                end
                                            end
                                        end
                                        init_meas_time=read_delay+1/input_vector_freq+2*in_vect_r_f_time;
                                        
                                        [MAX_settle_min_val,MAX_settle_min_idx]=max(settle_min(:));
                                        [MAX_settle_min_idx_part,MAX_settle_min_idx_col,MAX_settle_min_idx_pol]=ind2sub(size(settle_min),MAX_settle_min_idx);
                                        coordinates_MAX_settle_min=[MAX_settle_min_idx_part,MAX_settle_min_idx_col,MAX_settle_min_idx_pol];
                                        
                                        [MAX_settle_max_val,MAX_settle_max_idx]=max(settle_max(:));
                                        [MAX_settle_max_idx_part,MAX_settle_max_idx_col,MAX_settle_max_idx_pol]=ind2sub(size(settle_max),MAX_settle_max_idx);
                                        coordinates_MAX_settle_max=[MAX_settle_max_idx_part,MAX_settle_max_idx_col,MAX_settle_max_idx_pol];

                                        [MAX_settle_neuron_min_val,MAX_settle_neuron_min_idx]=max(settle_neuron_min(:));
                                        [MAX_settle_neuron_min_idx_part,MAX_settle_neuron_min_idx_col,MAX_settle_neuron_min_idx_pol]=ind2sub(size(settle_neuron_min),MAX_settle_neuron_min_idx);
                                        coordinates_MAX_settle_neuron_min=[MAX_settle_neuron_min_idx_part,MAX_settle_neuron_min_idx_col,MAX_settle_neuron_min_idx_pol];
                                        
                                        [MAX_settle_neuron_max_val,MAX_settle_neuron_max_idx]=max(settle_neuron_max(:));
                                        [MAX_settle_neuron_max_idx_part,MAX_settle_neuron_max_idx_col,MAX_settle_neuron_max_idx_pol]=ind2sub(size(settle_neuron_max),MAX_settle_neuron_max_idx);
                                        coordinates_MAX_settle_neuron_max=[MAX_settle_neuron_max_idx_part,MAX_settle_neuron_max_idx_col,MAX_settle_neuron_max_idx_pol];
                                        
                                        
                                        max_settling_time=max([MAX_settle_max_val MAX_settle_min_val])-init_meas_time;
                                        max_settling_neuron_time=max([MAX_settle_neuron_max_val MAX_settle_neuron_min_val])-init_meas_time;
                                        save(fullfile(results_folder,'RE_latency.mat'),'settle_max','settle_min','max_settling_time','settle_neuron_max','settle_neuron_min','max_settling_time','max_settling_neuron_time','init_meas_time','C_memdiode','C_interline', 'coordinates_MAX_settle_min','coordinates_MAX_settle_max','coordinates_MAX_settle_neuron_min','coordinates_MAX_settle_neuron_max');
                                    end
                                end
                            end
                        end
                    end
                else
                    raw_data_filename=strcat('NW_',num2str(neu_per_layer(1)),'by',num2str(neu_per_layer(2)),'_closed_loop.raw');
                    raw_data=fullfile(root_str,'netlists',raw_data_filename);
                    if exist(raw_data,'file')
                        movefile(raw_data,results_folder);
                    end

                    op_data_filename=strcat('NW_',num2str(neu_per_layer(1)),'by',num2str(neu_per_layer(2)),'_closed_loop.op.raw');
                    op_data=fullfile(root_str,'netlists',op_data_filename);
                    if exist(op_data,'file')
                        movefile(op_data,results_folder);
                    end

                    log_data_filename=strcat('NW_',num2str(neu_per_layer(1)),'by',num2str(neu_per_layer(2)),'_closed_loop.log');        
                    log_data=fullfile(root_str,'netlists',log_data_filename);
                    if exist(log_data,'file')
                        movefile(log_data,results_folder);
                    end
                    
                    RAW_DATA=LTspice2Matlab(fullfile(results_folder,raw_data_filename));
                    RAW_DATA.neu_per_layer=neu_per_layer;
                    RAW_DATA.simEngine='LTspice';
                    
                    save(fullfile(results_folder,'RAW_DATA.mat'),'RAW_DATA','W_matrix');

                    for i=1:size(RAW_DATA.variable_name_list,2)
                        variable_name=RAW_DATA.variable_name_list{i};
                        string_search=sprintf('V(write_ok!)');
                        if strcmp(string_search,variable_name)
                            WR_OK_pos=i;
                            break
                        end
                    end

                    if sim_WR==1
                        if RAW_DATA.num_steps>1
                            WR_OK_status=squeeze(RAW_DATA.variable_mat(WR_OK_pos,:,:));
                            diff_WR_OK_status=diff(WR_OK_status,1,1);
                            trans_WR_OK=find(diff_WR_OK_status>0);
                        else
                            WR_OK_status=RAW_DATA.variable_mat(WR_OK_pos,:);
                            diff_WR_OK_status=diff(WR_OK_status,1,2);
                            trans_WR_OK=find(diff_WR_OK_status>0);
                        end

                        if isempty(trans_WR_OK)==0
                            [idx_row,idx_col]=ind2sub(size(diff_WR_OK_status),trans_WR_OK);
                            if max(idx_col)<RAW_DATA.num_steps
                                simulation_time_aux=simulation_time*1.5;
                                SearchString=sprintf('.tran 0 %.12e 0 %.12e',simulation_time,max_time_step);
                                ReplaceString=sprintf('.tran 0 %.12e 0 %.12e',simulation_time_aux,max_time_step);
                                func_replace_string(fullfile(sim_folder,netlist_name), fullfile(sim_folder,netlist_name), SearchString, ReplaceString)
                                simulation_time=simulation_time_aux;
                            else
                                sim_complete=1;
                            end
                        else
                            simulation_time_aux=simulation_time*1.5;
                            SearchString=sprintf('.tran 0 %.12e 0 %.12e',simulation_time,max_time_step);
                            ReplaceString=sprintf('.tran 0 %.12e 0 %.12e',simulation_time_aux,max_time_step);
                            func_replace_string(fullfile(sim_folder,netlist_name), fullfile(sim_folder,netlist_name), SearchString, ReplaceString)
                            simulation_time=simulation_time_aux;                    
                        end
                    else
                        sim_complete=1;
                    end
                end
            end
            fprintf('done!\n');
            tElapsed_sim=toc(tStart_sim);
            
            if guiStatus_data
                str_to_display{end+1}=sprintf("%sTotal time required for the read/write simulation: %s",str_header,duration([0, 0, tElapsed_sim]));
                figure(status_gui);
                status_guiObject.log.String=str_to_display;
                map.cmd_history=str_to_display;
                status_guiObject.log.Value=length(str_to_display)+1;
                drawnow
            end
        end
        if strcmp(genPlot,'yes')
            if strcmp(analysis,'o_point')
            
            else
                if strcmpi(output_fmt,'ascii')
                    fprintf('---> Generating plots...')
                    processed_data=generate_plots_dualPart(RAW_DATA,partitions_N,1,1,1,0,1,results_folder,'networkSize',RAW_DATA.neu_per_layer);
                    save(fullfile(results_folder,'RAW_DATA.mat'),'processed_data','-append');
                end
            end
        end
        fprintf('done!\n');
end

