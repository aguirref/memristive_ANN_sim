function [data_sim] = fit_memdiode_lineal(model_name, read_voltage, varargin)
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
    
    if mod(nargin-2,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-2
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
        end
    end   

    file_Vt=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_data_set.txt');
    file_It=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_set.txt');
    file_It_2=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_set_2.txt');
    file_It_scaled=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_set_scaled.txt');
    file_Vt_reset=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_data_reset.txt');
    file_It_reset=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_reset.txt');
    file_It_reset_2=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_reset_2.txt');
    file_It_reset_scaled=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_reset_scaled.txt');   
    file_Vt_complete=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_complete.txt');
    
    hspice_root_dir='/usr/eda_tools/synopsys/HSPICE/O-2018.09/hspice/';

%     Vt_set=importdata(file_Vt);
%     It_set=importdata(file_It);
%     Vt_reset=importdata(file_Vt_reset);
%     It_reset=importdata(file_It_reset);   
%     It_reset(:,2)=It_reset(:,2)*(-1);
%     It_reset(:,1)=It_reset(:,1)+It_set(end,1);
%     Vt_reset(:,1)=Vt_reset(:,1)+Vt_set(end,1);
%     
%     Vt=vertcat(Vt_set,Vt_reset);
%     It=vertcat(It_set,It_reset);

    
    if strcmpi(region,'SET') || strcmpi(region,'both')
        Vt=importdata(file_Vt);
        It=importdata(file_It);
        if strcmpi(kR,'yes')
            It(:,2)=It(:,2).*scale;
        else
            if strcmpi(r2sweep,'ROFF')
                [max_val,max_pos]=max(Vt(:,2));
                It(1:floor((max_pos-1)/2)-1,2)=It(1:floor((max_pos-1)/2)-1,2).*scale;   
                scaling_vector=linspace(scale,1,(max_pos-1)-floor((max_pos-1)/2)+1)';
                It(floor((max_pos-1)/2):max_pos-1,2)=It(floor((max_pos-1)/2):max_pos-1,2).*scaling_vector;
            else
                [max_val,max_pos]=max(Vt(:,2));
                It(max_pos:end,2)=It(max_pos:end,2).*scale;
                It(1:floor((max_pos-1)/2)-1,2)=It(1:floor((max_pos-1)/2)-1,2).*scale_Roff;   
                scaling_vector=linspace(scale_Roff,scale,(max_pos-1)-floor((max_pos-1)/2)+1)';
                It(floor((max_pos-1)/2):max_pos-1,2)=It(floor((max_pos-1)/2):max_pos-1,2).*scaling_vector;                
            end
        end
    end    
    if strcmpi(region,'RESET') || strcmpi(region,'both')
        Vt_reset=importdata(file_Vt_reset);
        It_reset=importdata(file_It_reset);
        if strcmpi(kR,'yes')
            It_reset(:,2)=It_reset(:,2).*scale;
        else
            if strcmpi(r2sweep,'ROFF')
                [min_val,min_pos]=min(Vt_reset(:,2));
                It_reset(min_pos:end,2)=It_reset(min_pos:end,2).*scale;
                %It_reset(floor((min_pos-1)/2)-1:end,2)=It_reset(floor((min_pos-1)/2)-1:end,2).*scale.*(-1);   
                scaling_vector=flip(logspace(log10(scale),log10(1),(min_pos-1)-floor((min_pos-1)/2)+1)').*(1);
                It_reset(floor((min_pos-1)/2):min_pos-1,2)=It_reset(floor((min_pos-1)/2):min_pos-1,2).*scaling_vector;  
            else
                [min_val,min_pos]=min(Vt_reset(:,2));
                It_reset(min_pos:end,2)=It_reset(min_pos:end,2).*scale_Roff;
                %It_reset(1:min_pos,2)=It_reset(1:min_pos,2).*scale;
                It_reset(1:floor((min_pos-1)/2),2)=It_reset(1:floor((min_pos-1)/2),2).*scale;                
                scaling_vector=flip(logspace(log10(scale_Roff),log10(scale),(min_pos-1)-floor((min_pos-1)/2))').*(1);
                It_reset(floor((min_pos-1)/2)+1:min_pos-1,2)=It_reset(floor((min_pos-1)/2)+1:min_pos-1,2).*scaling_vector;                  
                
%                 scaling_vector=flip(logspace(log10(scale),log10(1),(min_pos-1)-floor((min_pos-1)/2)+1)').*(1);
%                 It_reset(floor((min_pos-1)/2):min_pos-1,2)=It_reset(floor((min_pos-1)/2):min_pos-1,2).*scaling_vector;  
%                 It(floor((max_pos-1)/2):max_pos-1,2)=It(floor((max_pos-1)/2):max_pos-1,2).*scaling_vector;   
            end
        end
    end
    
    
    if strcmpi(region,'SET')
        It=It;
        Vt=Vt; 
    elseif strcmpi(region,'RESET')
        It=It_reset;
        Vt=Vt_reset;           
    elseif strcmpi(region,'both')
        It_reset(:,2)=It_reset(:,2)*(-1);
        It_reset(:,1)=It_reset(:,1)+It(end,1);
        Vt_reset(:,1)=Vt_reset(:,1)+Vt(end,1);
        Vt=vertcat(Vt,Vt_reset);
        It=vertcat(It,It_reset);          
    end

    dlmwrite(file_It_scaled,It,'\t');
    dlmwrite(file_Vt_complete,Vt,'\t'); 
    
    RUNLVL=5;
    
    figure();
    ax1=subplot(1,2,1);
    semilogy(Vt(:,2),abs(It(:,2)),'k');
    ylabel('Memdiode current [A]');
    xlabel('Memdiode Voltage [V]');
    axis([-3 1.5 1e-7 1e-1]);
    hold on

    ax2=subplot(1,2,2);
    plot(Vt(:,2),abs(It(:,2)),'k');
    ylabel('Memdiode current [A]');
    xlabel('Memdiode Voltage [V]');
    axis([-3 1.5 0 1.2e-2]);
    hold on
    
    sim_V=Vt(:,2)';
    sim_I=It(:,2)';

    [~,aux_var_pos]=findpeaks(sim_V);

    sim_V_HRS=sim_V(1:aux_var_pos);
    sim_I_HRS=sim_I(1:aux_var_pos);
    aux_pos_minor=find(sim_V_HRS<read_voltage);
    aux_pos_major=find(sim_V_HRS>read_voltage);
    IMIN=interp1(sim_V_HRS([max(aux_pos_minor) min(aux_pos_major)]),sim_I_HRS([max(aux_pos_minor) min(aux_pos_major)]),read_voltage);

    sim_V_LRS=sim_V(aux_var_pos:end);
    sim_I_LRS=sim_I(aux_var_pos:end);
    aux_pos_minor=find(sim_V_LRS>read_voltage);
    aux_pos_major=find(sim_V_LRS<read_voltage);
    IMAX=interp1(sim_V_LRS([max(aux_pos_minor) min(aux_pos_major)]),sim_I_LRS([max(aux_pos_minor) min(aux_pos_major)]),read_voltage);
    
    ax1=subplot(1,2,1);
    semilogy(sim_V,sim_V*IMAX/read_voltage,'b');
    semilogy(sim_V,sim_V*IMIN/read_voltage,'r');
    
    %% Esta parte escribe el archivo del modelo para HSPICE
    
    if ~exist(fullfile(prj_dir,filesep,'models/memdiode/lineal'),'dir')
        mkdir(fullfile(prj_dir,filesep,'models/memdiode/lineal'))
    end
    model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/lineal/lineal_HSPICE_",model_name,'_',strrep(num2str(read_voltage),'.','p'),'.sp'));
    
    model_fileID = fopen(model_fileID_str,'w');
    fprintf(model_fileID,'***************\n');
    fprintf(model_fileID,'*For read_voltage=%f\n',read_voltage);
    
    fprintf(model_fileID,'.global gnd!\n');
    fprintf(model_fileID,'.subckt memdiode p n H0=0\n');
    fprintf(model_fileID,'.param gon=%.3e\n', IMAX/read_voltage);
    fprintf(model_fileID,'.param goff=%.3e\n', IMIN/read_voltage);

    fprintf(model_fileID,'*I-V\n');
    fprintf(model_fileID,' RS p n R=''1/(gon*H0+(1-H0)*goff)''\n');

    fprintf(model_fileID,'.ends memdiode\n');
    fclose(model_fileID);    
    
end