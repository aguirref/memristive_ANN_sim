function [data_sim] = fit_memdiode(results_folder, netlist_folder, version, model_name, varargin)
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
    gam=0;
    isb=1;
    H0_limit=2e-6;
    region='both';
    electroformed='no';
    r2sweep='RON';
    variability_model='no';
    scale_Roff=1e-3;
    mode_data='singleFile';
    
    hspice_root_dir='/usr/eda_tools/synopsys/HSPICE/O-2018.09/hspice/';
    
    prj_dir=fullfile('/home/users','aguirref','nn_rs_uab');
    root_folder=fullfile(prj_dir,'results','memdiode_sim');
    
    file_Vt=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_data_set.txt');
    file_It=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_set.txt');
    file_It_2=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_set_2.txt');
    file_It_scaled=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_set_scaled.txt');
    file_Vt_reset=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_data_reset.txt');
    file_It_reset=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_reset.txt');
    file_It_reset_2=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_reset_2.txt');
    file_It_reset_scaled=fullfile(prj_dir,'models','memdiode','model_fit','mean_It_data_reset_scaled.txt');   
    file_Vt_complete=fullfile(prj_dir,'models','memdiode','model_fit','mean_Vt_complete.txt');
    
    
    if mod(nargin-4,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-4
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
            if strcmp(varargin{i},'vtOpt')
                vtOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'VOrOpt')
                VOrOpt=varargin{i+1};
            end
            if strcmp(varargin{i},'VOsOpt')
                VOsOpt=varargin{i+1};
            end            
            if strcmp(varargin{i},'TOrOpt')
                TOrOpt=varargin{i+1};
            end                     
            if strcmp(varargin{i},'TOsOpt')
                TOsOpt=varargin{i+1};
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
            if strcmp(varargin{i},'H0_limit')
                H0_limit=varargin{i+1};
            end  
            if strcmp(varargin{i},'isb')
                isb=varargin{i+1};
            end  
            if strcmp(varargin{i},'gam')
                gam=varargin{i+1};
            end              
            if strcmp(varargin{i},'varModel')
                variability_model=varargin{i+1};
            end  
            if strcmp(varargin{i},'file_Vt')
                file_Vt=varargin{i+1};
            end              
            if strcmp(varargin{i},'file_It')
                file_It=varargin{i+1};
            end  
        end
    end   

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

if strcmpi(mode_data,'singleFile')
    Vt=importdata(file_Vt);
    It=importdata(file_It);
    if strcmpi(kR,'yes')
        It(:,2)=It(:,2).*scale;
    end
else
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
    
    netlist_name='memdiode_fit_netlist.sp';
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
    
    if strcmpi(version,'new')
        fprintf(fileID,'.inc ''%s/models/memdiode/model_fit/memdiode_HSPICE_new_base.sp''\n',prj_dir);
        fprintf(fileID,'Xmemdiodo pos gnd! h memdiode H0=''H_init'' beta=''beta_opt'' T0s=''T0s_opt'' V0s=''V0s_opt'' T0r=''T0r_opt'' V0r=''V0r_opt'' imax=''imax_opt'' alphamax=''alphamax_opt'' rsmax=''rs_opt'' imin=''imin_opt'' alphamin=''alphamin_opt'' rsmin=''rs_opt''\n');
    elseif strcmpi(version,'old_SBSF')
        fprintf(fileID,'.inc ''%s/models/memdiode/model_fit/memdiode_HSPICE_old_SBSF_base.sp''\n',prj_dir);
        fprintf(fileID,'Xmemdiodo pos gnd! h memdiode H0=''H_init'' H0_limit=%.3e beta=''beta_opt'' etaset=''etaset_opt'' vset=''vset_opt'' vt=''vt_opt'' etares=''etares_opt'' vres=''vres_opt'' imax=''imax_opt'' alphamax=''alphamax_opt'' rsmax=''rsmax_opt'' imin=''imin_opt'' alphamin=''alphamin_opt'' rsmin=''rsmin_opt'' isb=''isb_opt'' gam=''gam_opt''\n', H0_limit);        
    else
        fprintf(fileID,'.inc ''%s/models/memdiode/model_fit/memdiode_HSPICE_old_base.sp''\n',prj_dir);
        fprintf(fileID,'Xmemdiodo pos gnd! h memdiode H0=''H_init'' beta=''beta_opt'' etaset=''etaset_opt'' vset=''vset_opt'' etares=''etares_opt'' vres=''vres_opt'' imax=''imax_opt'' alphamax=''alphamax_opt'' rsmax=''rsmax_opt'' imin=''imin_opt'' alphamin=''alphamin_opt'' rsmin=''rsmin_opt''\n');        
    end
    fprintf(fileID,'Rh h gnd! R=100k\n');
    fprintf(fileID,'Vsignal pos gnd! PWL PWLFILE=''%s''\n',file_Vt_complete);
    %fprintf(fileID,'Rcomp   pos-1 pos  R=''I(Rcomp)<1e-3 ? 0.1 : V(pos-1,pos)/1e-3''\n');
    
    if strcmpi(version,'new')
        fprintf(fileID,'.param\n');
        fprintf(fileID,'+ H_init            = 0\n');
        fprintf(fileID,'+ beta_opt    		= 0.5\n');%= opt1(0.5, 0.4, 0.6)\n');
        fprintf(fileID,'+ T0s_opt     		= opt1(%.6e,%.6e,%.6e)\n',TOsOpt(1),TOsOpt(2),TOsOpt(3));
        fprintf(fileID,'+ V0s_opt    		= opt1(%.6e,%.6e,%.6e)\n',VOsOpt(1),VOsOpt(2),VOsOpt(3));
        fprintf(fileID,'+ T0r_opt     		= opt1(%.6e,%.6e,%.6e)\n',TOrOpt(1),TOrOpt(2),TOrOpt(3));
        fprintf(fileID,'+ V0r_opt     		= opt1(%.6e,%.6e,%.6e)\n',VOrOpt(1),VOrOpt(2),VOrOpt(3));
        fprintf(fileID,'+ imax_opt     		= opt1(%.6e,%.6e,%.6e)\n',imaxOpt(1),imaxOpt(2),imaxOpt(3));
        fprintf(fileID,'+ alphamax_opt    	= opt1(%.6e,%.6e,%.6e)\n',alphamaxOpt(1),alphamaxOpt(2),alphamaxOpt(3));
        fprintf(fileID,'+ rsmax_opt      	= rs_opt\n');%= opt1(50,1,200)\n');
        fprintf(fileID,'+ imin_opt   		= opt1(%.6e,%.6e,%.6e)\n',iminOpt(1),iminOpt(2),iminOpt(3));
        fprintf(fileID,'+ alphamin_opt      = opt1(%.6e,%.6e,%.6e)\n',alphaminOpt(1),alphaminOpt(2),alphaminOpt(3));
        fprintf(fileID,'+ rsmin_opt   		= rs_opt\n');%= opt1(50,1,200)\n\n');
        fprintf(fileID,'+ rs_opt   		    = opt1(%.6e,%.6e,%.6e)\n',rsOpt(1),rsOpt(2),rsOpt(3));
    else
        fprintf(fileID,'.param\n');
        if strcmpi(version,'old')
            fprintf(fileID,'+ H_init            = 0\n');
        elseif strcmpi(version,'old_SBSF')
            fprintf(fileID,'+ H_init            = 1.0e-5\n');
            if length(vtOpt)>1
                fprintf(fileID,'+ vt_opt    		= opt1(%.6e,%.6e,%.6e)\n',vtOpt(1),vtOpt(2),vtOpt(3));
            else
                fprintf(fileID,'+ vt_opt    		= %.6e\n',vtOpt);
            end
            fprintf(fileID,'+ gam_opt    		= %.6e\n',gam);
            fprintf(fileID,'+ isb_opt    		= %.6e\n',isb);
        end
        fprintf(fileID,'+ beta_opt    		= 0.5\n');
        fprintf(fileID,'+ etaset_opt   		= opt1(%.6e,%.6e,%.6e)\n',etasetOpt(1),etasetOpt(2),etasetOpt(3));
        if length(vsetOpt)>1
            fprintf(fileID,'+ vset_opt    		= opt1(%.6e,%.6e,%.6e)\n',vsetOpt(1),vsetOpt(2),vsetOpt(3));
        else
            fprintf(fileID,'+ vset_opt    		= %.6e\n',vsetOpt);
        end
        fprintf(fileID,'+ etares_opt     	= opt1(%.6e,%.6e,%.6e)\n',etaresOpt(1),etaresOpt(2),etaresOpt(3));
        fprintf(fileID,'+ vres_opt     		= opt1(%.6e,%.6e,%.6e)\n',vresOpt(1),vresOpt(2),vresOpt(3));
        fprintf(fileID,'+ imax_opt     		= opt1(%.6e,%.6e,%.6e)\n',imaxOpt(1),imaxOpt(2),imaxOpt(3));
        fprintf(fileID,'+ alphamax_opt    	= opt1(%.6e,%.6e,%.6e)\n',alphamaxOpt(1),alphamaxOpt(2),alphamaxOpt(3));
        fprintf(fileID,'+ rsmax_opt      	= rs_opt\n');%opt1(%.6e,%.6e,%.6e)\n',rsmaxOpt(1),rsmaxOpt(2),rsmaxOpt(3));   
        fprintf(fileID,'+ imin_opt   		= opt1(%.6e,%.6e,%.6e)\n',iminOpt(1),iminOpt(2),iminOpt(3));
        fprintf(fileID,'+ alphamin_opt      = opt1(%.6e,%.6e,%.6e)\n',alphaminOpt(1),alphaminOpt(2),alphaminOpt(3));
        fprintf(fileID,'+ rsmin_opt   		= rs_opt\n');%= opt1(%.6e,%.6e,%.6e)\n',rsminOpt(1),rsminOpt(2),rsminOpt(3));   
        fprintf(fileID,'+ CH0_opt   		= %.6e\n',CH0_opt);%= opt1(50,1,200)\n\n');
        fprintf(fileID,'+ rs_opt   		    = opt1(%.6e,%.6e,%.6e)\n',rsOpt(1),rsOpt(2),rsOpt(3));        
    end

    fprintf(fileID,'*data\n');
    fprintf(fileID,'.param imem=0\n');
    fprintf(fileID,'.data memdata MER\n');
    fprintf(fileID,'+ FILE=''%s''\n',file_It_scaled);
    fprintf(fileID,'+ TIME=1 imem=2\n');
    fprintf(fileID,'+ OUT=''%s''\n',file_It_2);
    fprintf(fileID,'.enddata\n');

    fprintf(fileID,'.PRINT TRAN V(pos) I(xmemdiodo.RS) V(h)\n');
    
    fprintf(fileID,'.tran 10m %f\n',Vt(end,1));

    fprintf(fileID,'*..optimization sweeps\n');
    fprintf(fileID,'.meas tran comp1 err2 par(imem) I(xmemdiodo.RS) minval=1e-15 ignor=1e-18\n');
    fprintf(fileID,'.tran DATA = memdata  SWEEP OPTIMIZE = opt1  RESULTS = comp1  MODEL = optmod\n');

    fprintf(fileID,'.end\n');

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
        disp(tline)
        if strcmp(tline, 'x'),'k'
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
    
    data_sim_pre=importdata(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0-pre')));
    data_sim_post=importdata(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0-post')));
    
    fid = fopen(fullfile(results_folder,strrep(netlist_name,'.sp','.lis')));
    %fid2 = fopen(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0-pre')),'w+');
    %fid3 = fopen(fullfile(results_folder,strrep(netlist_name,'.sp','.printtr0-post')),'w+');
    
    tline = fgetl(fid);
    i_param = 1;
    while ischar(tline)
        disp(tline)
        if strfind(tline, '.param beta_opt')
            A=split(tline);
            params{i_param,1}=A{3,1};
            params{i_param,2}=A{5,1};
            i_param=i_param+1;
            beta_opt=A{5,1};
        end
        if strcmpi(version,'new')
            if strfind(tline, '.param t0s_opt')
                A=split(tline);
                params{i_param,1}=A{3,1};
                params{i_param,2}=A{5,1};
                i_param=i_param+1;
                t0s_opt=A{5,1};
            end
            if strfind(tline, '.param v0s_opt')
                A=split(tline);
                params{i_param,1}=A{3,1};
                params{i_param,2}=A{5,1};
                i_param=i_param+1;
                v0s_opt=A{5,1};
            end
            if strfind(tline, '.param t0r_opt')
                A=split(tline);
                params{i_param,1}=A{3,1};
                params{i_param,2}=A{5,1};
                i_param=i_param+1;
                t0r_opt=A{5,1};
            end
            if strfind(tline, '.param v0r_opt')
                A=split(tline);
                params{i_param,1}=A{3,1};
                params{i_param,2}=A{5,1};
                i_param=i_param+1;
                vOr_opt=A{5,1};
            end
        else
            if strfind(tline, '.param etaset_opt')
                A=split(tline);
                params{i_param,1}=A{3,1};
                params{i_param,2}=A{5,1};
                i_param=i_param+1;
                etaset_opt=A{5,1};
            end
            if strfind(tline, '.param vset_opt')
                A=split(tline);
                params{i_param,1}=A{3,1};
                params{i_param,2}=A{5,1};
                i_param=i_param+1;
                vset_opt=A{5,1};
            end
            if strfind(tline, '.param vt_opt')
                A=split(tline);
                params{i_param,1}=A{3,1};
                params{i_param,2}=A{5,1};
                i_param=i_param+1;
                vt_opt=A{5,1};
            end
            if strfind(tline, '.param etares_opt')
                A=split(tline);
                params{i_param,1}=A{3,1};
                params{i_param,2}=A{5,1};
                i_param=i_param+1;
                etares_opt=A{5,1};
            end
            if strfind(tline, '.param vres_opt')
                A=split(tline);
                params{i_param,1}=A{3,1};
                params{i_param,2}=A{5,1};
                i_param=i_param+1;
                vres_opt=A{5,1};
            end            
        end
        if strfind(tline, '.param imax_opt')
            A=split(tline);
            params{i_param,1}=A{3,1};
            params{i_param,2}=A{5,1};
            i_param=i_param+1;
            imax_opt=A{5,1};
        end
        if strfind(tline, '.param alphamax_opt')
            A=split(tline);
            params{i_param,1}=A{3,1};
            params{i_param,2}=A{5,1};
            i_param=i_param+1;
            alphamax_opt=A{5,1};
        end        
        if strfind(tline, '.param rsmax_opt')
            A=split(tline);
            params{i_param,1}=A{3,1};
            params{i_param,2}=A{5,1};
            i_param=i_param+1;
            rsmax_opt=A{5,1};
        end        
        if strfind(tline, '.param imin_opt')
            A=split(tline);
            params{i_param,1}=A{3,1};
            params{i_param,2}=A{5,1};
            i_param=i_param+1;
            imin_opt=A{5,1};
        end  
        if strfind(tline, '.param alphamin_opt')
            A=split(tline);
            params{i_param,1}=A{3,1};
            params{i_param,2}=A{5,1};
            i_param=i_param+1;
            alphamin_opt=A{5,1};
        end        
        if strfind(tline, '.param rsmin_opt')
            A=split(tline);
            params{i_param,1}=A{3,1};
            params{i_param,2}=A{5,1};
            i_param=i_param+1;
            rsmin_opt=A{5,1};
        end   
        if strfind(tline, '.param rs_opt')
            A=split(tline);
            params{i_param,1}=A{3,1};
            params{i_param,2}=A{5,1};
            i_param=i_param+1;
            rs_opt=A{5,1};
        end   
        if strfind(tline, '.param CH0_opt')
            A=split(tline);
            params{i_param,1}=A{3,1};
            params{i_param,2}=A{5,1};
            i_param=i_param+1;
            CH0_opt=A{5,1};
        end           
        tline = fgetl(fid);
        lineCounter = lineCounter + 1;                
    end

    data_sim.data_fit=params;
    data_sim.data_pre_fit=data_sim_pre;
    data_sim.data_post_fit=data_sim_post;
    
    sim_I_pre=data_sim_pre.data(:,3);
    sim_time_pre=data_sim_pre.data(:,1);
    sim_V_pre=data_sim_pre.data(:,2);
    sim_h_pre=data_sim_pre.data(:,4);
    
    
    sim_I_post=data_sim_post.data(:,3);
    sim_time_post=data_sim_post.data(:,1);
    sim_V_post=data_sim_post.data(:,2);
    sim_h_post=data_sim_post.data(:,4);
    
    axes(ax1);
    semilogy(sim_V_pre,abs(sim_I_pre),'--r');
    semilogy(sim_V_post,abs(sim_I_post),'-r');
    semilogy(sim_V_pre,abs(sim_h_pre),'--b');
    semilogy(sim_V_post,abs(sim_h_post),'-b');
    
    axes(ax2);
    plot(sim_V_pre,abs(sim_I_pre),'--r');    
    plot(sim_V_post,abs(sim_I_post),'-r'); 
    plot(sim_V_pre,abs(sim_h_pre),'--b');    
    plot(sim_V_post,abs(sim_h_post),'-b'); 
    %% Esta parte escribe el archivo del modelo para HSPICE
    if strcmpi(version,'new')
        if strcmpi(variability_model,'yes')
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/new_version/memdiode_HSPICE_new_",model_name,'_var.sp'));
        else
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/new_version/memdiode_HSPICE_new_",model_name,'.sp'));
        end
    elseif strcmpi(version,'old_SBSF')
        if strcmpi(variability_model,'yes')
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/old_SBSF_version/memdiode_HSPICE_old_SBSF_",model_name,'_var.sp'));
        else
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/old_SBSF_version/memdiode_HSPICE_old_SBSF_",model_name,'.sp'));
        end
    else
        if strcmpi(variability_model,'yes')
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/old_version/memdiode_HSPICE_old_",model_name,'_var.sp'));
        else
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/old_version/memdiode_HSPICE_old_",model_name,'.sp'));
        end
    end
    model_fileID = fopen(model_fileID_str,'w');
    if strcmpi(version,'new')
        fprintf(model_fileID,'***************\n');
        fprintf(model_fileID,'.global gnd!\n');
        fprintf(model_fileID,'.subckt memdiode p n H0=0\n');
        fprintf(model_fileID,'.param beta=%.12e\n',beta_opt);
        fprintf(model_fileID,'.param EI=1e-12\n');
        fprintf(model_fileID,'.param T0s=%.12e\n',str2double(t0s_opt));
        fprintf(model_fileID,'.param V0s=%.12e\n',str2double(v0s_opt));
        fprintf(model_fileID,'.param T0r=%.12e\n',t0r_opt);
        fprintf(model_fileID,'.param V0r=%.12e\n',v0r_opt);
        fprintf(model_fileID,'.param imax=%.12e\n',str2double(imax_opt));
        fprintf(model_fileID,'.param alphamax=%.12e\n',str2double(alphamax_opt));
        fprintf(model_fileID,'.param rsmax=%.12e\n',str2double(rs_opt));
        fprintf(model_fileID,'.param imin=%.12e\n',str2double(imin_opt));
        fprintf(model_fileID,'.param alphamin=%.12e\n',str2double(alphamin_opt));
        fprintf(model_fileID,'.param rsmin=%.12e\n\n',str2double(rs_opt));

        fprintf(model_fileID,'*Auxiliary functions\n');
        fprintf(model_fileID,'.param I0(x)=''imax*x+imin*(1-x)''\n');
        fprintf(model_fileID,'.param A(x)=''alphamax*x+alphamin*(1-x)''\n');
        fprintf(model_fileID,'.param Rss(x)=''rsmax*x+rsmin*(1-x)''\n\n');

        fprintf(model_fileID,'*H-V\n');
        fprintf(model_fileID,'EV A gnd! vol=1\n');
        fprintf(model_fileID,'RH H A    R=''T0s*exp(-V(p,n)/V0s)''\n');
        fprintf(model_fileID,'RD H gnd! R=''T0r*exp(V(p,n)/V0r)''\n');
        fprintf(model_fileID,'CH H gnd! C=1 IC=''H0''\n\n');

        fprintf(model_fileID,'*I-V\n');
        fprintf(model_fileID,'RS p D R=''Rss(V(H))''\n');
        fprintf(model_fileID,'GD D n cur=''I0(V(H))*(exp(beta*A(V(H))*V(D,n))-exp(-(1-beta)*A(V(H))*V(D,n)))+EI''\n\n');

        fprintf(model_fileID,'.ends memdiode\n');
    elseif strcmpi(version,'old')
        fprintf(model_fileID,'***************\n');
        fprintf(model_fileID,'.global gnd!\n');
        if strcmpi(variability_model,'yes')
            fprintf(model_fileID,'.subckt memdiode p n H0=0\n');
            fprintf(model_fileID,'+ imax=%.12e\n',str2double(imax_opt));
            fprintf(model_fileID,'+ alphamax=%.12e\n',str2double(alphamax_opt));
            fprintf(model_fileID,'+ rsmax=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'+ imin=%.12e\n',str2double(imin_opt));
            fprintf(model_fileID,'+ alphamin=%.12e\n',str2double(alphamin_opt));
            fprintf(model_fileID,'+ rsmin=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'+ etaset=%.12e\n',str2double(etaset_opt));
            fprintf(model_fileID,'+ vset=%.12e\n',str2double(vset_opt));
            fprintf(model_fileID,'+ etares=%.12e\n',str2double(etares_opt));
            fprintf(model_fileID,'+ vres=%.12e\n\n',str2double(vres_opt));
        else
            fprintf(model_fileID,'.subckt memdiode p n H0=0\n');
            fprintf(model_fileID,'.param imax=%.12e\n',str2double(imax_opt));
            fprintf(model_fileID,'.param alphamax=%.12e\n',str2double(alphamax_opt));
            fprintf(model_fileID,'.param rsmax=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'.param imin=%.12e\n',str2double(imin_opt));
            fprintf(model_fileID,'.param alphamin=%.12e\n',str2double(alphamin_opt));
            fprintf(model_fileID,'.param rsmin=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'.param etaset=%.12e\n',str2double(etaset_opt));
            fprintf(model_fileID,'.param vset=%.12e\n',str2double(vset_opt));
            fprintf(model_fileID,'.param etares=%.12e\n',str2double(etares_opt));
            fprintf(model_fileID,'.param vres=%.12e\n',str2double(vres_opt));
        end
        fprintf(model_fileID,'.param CH0=%.12e\n',CH0_opt);
        fprintf(model_fileID,'.param beta=%.12e\n',beta_opt);
        fprintf(model_fileID,'.param EI=0\n\n');

        fprintf(model_fileID,'*Auxiliary functions\n');
        fprintf(model_fileID,'.param I0(x)=''imax*x+imin*(1-x)''\n');
        fprintf(model_fileID,'.param A(x)=''alphamax*x+alphamin*(1-x)''\n');
        fprintf(model_fileID,'.param RS(x)=''rsmax*x+rsmin*(1-x)''\n');
        fprintf(model_fileID,'.param R(x)=''1/(1+exp(-etares*(x-vres)))''\n');
        fprintf(model_fileID,'.param S(x)=''1/(1+exp(-etaset*(x-vset)))''\n\n');

        fprintf(model_fileID,'*H-V\n');
        fprintf(model_fileID,'GH gnd! H cur=''min(R(V(p,n)),max(S(V(p,n)),V(H)))''\n');
        fprintf(model_fileID,'Rpar gnd! H R=1\n');
        fprintf(model_fileID,'CH H gnd! C=''CH0'' IC=''H0''\n\n');

        fprintf(model_fileID,'*I-V\n');
        fprintf(model_fileID,'RS p D R=''RS(V(H))''\n');
%        fprintf(model_fileID,'GD D n cur=''sgn(V(D,n))*I0(V(H))*(exp(beta*A(V(H))*abs(V(D,n)))-exp(-(1-beta)*A(V(H))*abs(V(D,n))))+EI''\n\n');
        fprintf(model_fileID,'GD D n cur=''I0(V(H))*(exp(beta*A(V(H))*V(D,n))-exp(-(1-beta)*A(V(H))*V(D,n)))+EI''\n\n');

        fprintf(model_fileID,'.ends memdiode\n'); 
    elseif strcmpi(version,'old_SBSF')
        fprintf(model_fileID,'***************\n');
        fprintf(model_fileID,'.global gnd!\n');
        if strcmpi(variability_model,'yes')
            fprintf(model_fileID,'.subckt memdiode p n H0=0\n');
            fprintf(model_fileID,'+ imax=%.12e\n',str2double(imax_opt));
            fprintf(model_fileID,'+ alphamax=%.12e\n',str2double(alphamax_opt));
            fprintf(model_fileID,'+ rsmax=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'+ imin=%.12e\n',str2double(imin_opt));
            fprintf(model_fileID,'+ alphamin=%.12e\n',str2double(alphamin_opt));
            fprintf(model_fileID,'+ rsmin=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'+ etaset=%.12e\n',str2double(etaset_opt));
            fprintf(model_fileID,'+ vset=%.12e\n',str2double(vset_opt));
            fprintf(model_fileID,'+ etares=%.12e\n',str2double(etares_opt));
            fprintf(model_fileID,'+ vres=%.12e\n\n',str2double(vres_opt));
        else
            fprintf(model_fileID,'.subckt memdiode p n H0=0\n');
            fprintf(model_fileID,'.param imax=%.12e\n',str2double(imax_opt));
            fprintf(model_fileID,'.param alphamax=%.12e\n',str2double(alphamax_opt));
            fprintf(model_fileID,'.param rsmax=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'.param imin=%.12e\n',str2double(imin_opt));
            fprintf(model_fileID,'.param alphamin=%.12e\n',str2double(alphamin_opt));
            fprintf(model_fileID,'.param rsmin=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'.param etaset=%.12e\n',str2double(etaset_opt));
            if length(vsetOpt)>1
                fprintf(model_fileID,'.param vset=%.12e\n',str2double(vset_opt));
            else
                fprintf(model_fileID,'.param vset=%.12e\n',vsetOpt);                
            end
            if length(vtOpt)>1
                fprintf(model_fileID,'.param vt=%.12e\n',str2double(vt_opt));
            else
                fprintf(model_fileID,'.param vt=%.12e\n',vtOpt);                
            end
            fprintf(model_fileID,'.param etares=%.12e\n',str2double(etares_opt));
            fprintf(model_fileID,'.param vres=%.12e\n',str2double(vres_opt));
        end
        fprintf(model_fileID,'.param CH0=%.12e\n',CH0_opt);
        fprintf(model_fileID,'.param beta=%.12e\n',beta_opt);
        fprintf(model_fileID,'.param EI=0\n\n');
        fprintf(model_fileID,'.param gam=%.3e\n', gam);
        fprintf(model_fileID,'.param gam0=0\n');
        fprintf(model_fileID,'.param isb=%.3e\n', isb);
        fprintf(model_fileID,'.param RPP=1e9\n');
        fprintf(model_fileID,'.param ri=1\n\n');

        fprintf(model_fileID,'*Auxiliary functions\n');
        fprintf(model_fileID,'.param I0(x)=''x>%.2e ? imin+(imax-imin)*x : 0''\n', H0_limit);
        fprintf(model_fileID,'.param A(x)=''alphamin+(alphamax-alphamin)*x''\n');
        fprintf(model_fileID,'.param RS(x)=''rsmin+(rsmax-rsmin)*x''\n');
        fprintf(model_fileID,'.param VSB(x)=''x>isb ? vt : vset''\n');
        fprintf(model_fileID,'.param ISF(x)=''gam==0 ? 1 : (pwr(x,gam)-gam0)''\n');
        fprintf(model_fileID,'.param S(x)=''1/(1+exp(-etaset*(x-VSB(I(RS)))))''\n');
        fprintf(model_fileID,'.param R(x)=''1/(1+exp(-etares*ISF(V(H))*(x-vres)))''\n');
        %fprintf(model_fileID,'.param limit(x)=''x<1 ? (x>0 ? x : 0 ) : 1''\n\n');
        
        fprintf(model_fileID,'*Quasi-static model: H-V\n');
        fprintf(model_fileID,'GH   gnd! H 	cur=''min(R(V(C,n)),max(S(V(C,n)),V(H)))''\n');
        fprintf(model_fileID,'Rpar gnd! H 	R=1\n');
        fprintf(model_fileID,'CH   H    gnd! 	C=''CH0'' IC=''H0''\n\n');

        fprintf(model_fileID,'*Quasi-static model: I-V\n');
        fprintf(model_fileID,'RE   p    C 	R=''ri''\n');
        fprintf(model_fileID,'RS   C 	B 	R=''RS(V(H))''\n');
        fprintf(model_fileID,'GD   B    n 	cur=''I0(V(H))*(exp(beta*A(V(H))*V(B,n))-exp(-(1-beta)*A(V(H))*V(B,n)))''\n');
        fprintf(model_fileID,'RB   p    n 	R=''RPP''\n');
        fprintf(model_fileID,'*GR   p    n 	cur=1e-12\n');     

        fprintf(model_fileID,'.ends memdiode\n');         
    end
    fclose(model_fileID);

    %% Esta parte escribe el archivo del modelo para LTSpice
    if strcmpi(version,'new')
        if strcmpi(variability_model,'yes')
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/new_version/memdiode_LTspice_new_",model_name,'_var.sp'));
        else
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/new_version/memdiode_LTspice_new_",model_name,'.sp'));
        end
    else
        if strcmpi(variability_model,'yes')
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/old_version/memdiode_LTspice_old_",model_name,'_var.sp'));
        else
            model_fileID_str=fullfile(prj_dir,filesep,strcat("models/memdiode/old_version/memdiode_LTspice_old_",model_name,'.sp'));
        end
    end
    model_fileID = fopen(model_fileID_str,'w');
    if strcmpi(version,'new')
        fprintf(model_fileID,'***************\n');
        fprintf(model_fileID,'.subckt memdiode p n H0=0\n');
        fprintf(model_fileID,'.param beta=%.12e\n',beta_opt);
        fprintf(model_fileID,'.param EI=0\n');
        fprintf(model_fileID,'.param T0s=%.12e\n',str2double(t0s_opt));
        fprintf(model_fileID,'.param V0s=%.12e\n',str2double(v0s_opt));
        fprintf(model_fileID,'.param T0r=%.12e\n',str2double(t0r_opt));
        fprintf(model_fileID,'.param V0r=%.12e\n',str2double(v0r_opt));
        fprintf(model_fileID,'.param imax=%.12e\n',str2double(imax_opt));
        fprintf(model_fileID,'.param alphamax=%.12e\n',str2double(alphamax_opt));
        fprintf(model_fileID,'.param rsmax=%.12e\n',str2double(rs_opt));
        fprintf(model_fileID,'.param imin=%.12e\n',str2double(imin_opt));
        fprintf(model_fileID,'.param alphamin=%.12e\n',str2double(alphamin_opt));
        fprintf(model_fileID,'.param rsmin=%.12e\n',str2double(rs_opt));

        fprintf(model_fileID,'*Auxiliary functions\n');
        fprintf(model_fileID,'.param I0(x)=''imax*x+imin*(1-x)''\n');
        fprintf(model_fileID,'.param A(x)=''alphamax*x+alphamin*(1-x)''\n');
        fprintf(model_fileID,'.param Rss(x)=''rsmax*x+rsmin*(1-x)''\n\n');

        fprintf(model_fileID,'*H-V\n');
        fprintf(model_fileID,'EV A 0 vol=1\n');
        fprintf(model_fileID,'RH H A R=''T0s*exp(-V(p,n)/V0s)''\n');
        fprintf(model_fileID,'RD H 0 R=''T0r*exp(V(p,n)/V0r)''\n');
        fprintf(model_fileID,'CH H 0 1 IC=''H0''\n\n');

        fprintf(model_fileID,'*I-V\n');
        fprintf(model_fileID,' RS p D R=''Rss(V(H))''\n');
        fprintf(model_fileID,'GD D n cur=''I0(V(H))*(exp(beta*A(V(H))*V(D,n))-exp(-(1-beta)*A(V(H))*V(D,n)))+EI''\n\n');

        fprintf(model_fileID,'.ends memdiode\n');
    elseif strcmpi(version,'old')
        fprintf(model_fileID,'***************\n');
        fprintf(model_fileID,'.global gnd!\n');
        if strcmpi(variability_model,'yes')
            fprintf(model_fileID,'.subckt memdiode + - H H0=0\n');
            fprintf(model_fileID,'+ imax=%.12e\n',str2double(imax_opt));
            fprintf(model_fileID,'+ alphamax=%.12e\n',str2double(alphamax_opt));
            fprintf(model_fileID,'+ rsmax=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'+ imin=%.12e\n',str2double(imin_opt));
            fprintf(model_fileID,'+ alphamin=%.12e\n',str2double(alphamin_opt));
            fprintf(model_fileID,'+ rsmin=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'+ etaset=%.12e\n',str2double(etaset_opt));
            fprintf(model_fileID,'+ vset=%.12e\n',str2double(vset_opt));
            fprintf(model_fileID,'+ etares=%.12e\n',str2double(etares_opt));
            fprintf(model_fileID,'+ vres=%.12e\n\n',str2double(vres_opt));
        else
            fprintf(model_fileID,'.subckt memdiode + - H H0=0\n');
            fprintf(model_fileID,'.param imax=%.12e\n',str2double(imax_opt));
            fprintf(model_fileID,'.param alphamax=%.12e\n',str2double(alphamax_opt));
            fprintf(model_fileID,'.param rsmax=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'.param imin=%.12e\n',str2double(imin_opt));
            fprintf(model_fileID,'.param alphamin=%.12e\n',str2double(alphamin_opt));
            fprintf(model_fileID,'.param rsmin=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'.param etaset=%.12e\n',str2double(etaset_opt));
            fprintf(model_fileID,'.param vset=%.12e\n',str2double(vset_opt));
            fprintf(model_fileID,'.param etares=%.12e\n',str2double(etares_opt));
            fprintf(model_fileID,'.param vres=%.12e\n',str2double(vres_opt));
        end
        fprintf(model_fileID,'.param CH0=%.12e\n',CH0_opt);
        fprintf(model_fileID,'.param beta=%.12e\n',beta_opt);
        fprintf(model_fileID,'.param EI=0\n\n');

        fprintf(model_fileID,'*Auxiliary functions\n');
        fprintf(model_fileID,'.func I0(x)=imax*x+imin*(1-x)\n');
        fprintf(model_fileID,'.func A(x)=alphamax*x+alphamin*(1-x)\n');
        fprintf(model_fileID,'.func RS(x)=rsmax*x+rsmin*(1-x)\n');
        fprintf(model_fileID,'.func R(x)=1/(1+exp(-etares*(x-vres)))\n');
        fprintf(model_fileID,'.func S(x)=1/(1+exp(-etaset*(x-vset)))\n\n');

        fprintf(model_fileID,'*H-V\n');
        fprintf(model_fileID,'BH 0 H I=min(R(V(+,-)),max(S(V(+,-)),V(H))) Rpar=1\n');
        fprintf(model_fileID,'CH H 0 {CH0} ic={H0}\n\n');
        
        fprintf(model_fileID,'*I-V\n');
        fprintf(model_fileID,'RS + D R=RS(V(H))\n');
%        fprintf(model_fileID,'GD D n cur=''sgn(V(D,n))*I0(V(H))*(exp(beta*A(V(H))*abs(V(D,n)))-exp(-(1-beta)*A(V(H))*abs(V(D,n))))+EI''\n\n');
        fprintf(model_fileID,'BD D - I=I0(V(H))*(exp(beta*A(V(H))*V(D,-))-exp(-(1-beta)*A(V(H))*V(D,-)))+EI\n\n');

        fprintf(model_fileID,'.ends memdiode\n'); 
    elseif strcmpi(version,'old_SBSF')
        fprintf(model_fileID,'***************\n');
        fprintf(model_fileID,'.global gnd!\n');
        if strcmpi(variability_model,'yes')
            fprintf(model_fileID,'.subckt memdiode p n H0=0\n');
            fprintf(model_fileID,'+ imax=%.12e\n',str2double(imax_opt));
            fprintf(model_fileID,'+ alphamax=%.12e\n',str2double(alphamax_opt));
            fprintf(model_fileID,'+ rsmax=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'+ imin=%.12e\n',str2double(imin_opt));
            fprintf(model_fileID,'+ alphamin=%.12e\n',str2double(alphamin_opt));
            fprintf(model_fileID,'+ rsmin=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'+ etaset=%.12e\n',str2double(etaset_opt));
            fprintf(model_fileID,'+ vset=%.12e\n',str2double(vset_opt));
            fprintf(model_fileID,'+ etares=%.12e\n',str2double(etares_opt));
            fprintf(model_fileID,'+ vres=%.12e\n\n',str2double(vres_opt));
        else
            fprintf(model_fileID,'.subckt memdiode p n H0=0\n');
            fprintf(model_fileID,'.param imax=%.12e\n',str2double(imax_opt));
            fprintf(model_fileID,'.param alphamax=%.12e\n',str2double(alphamax_opt));
            fprintf(model_fileID,'.param rsmax=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'.param imin=%.12e\n',str2double(imin_opt));
            fprintf(model_fileID,'.param alphamin=%.12e\n',str2double(alphamin_opt));
            fprintf(model_fileID,'.param rsmin=%.12e\n',str2double(rs_opt));
            fprintf(model_fileID,'.param etaset=%.12e\n',str2double(etaset_opt));
            if length(vsetOpt)>1
                fprintf(model_fileID,'.param vset=%.12e\n',str2double(vset_opt));
            else
                fprintf(model_fileID,'.param vset=%.12e\n',vsetOpt);                
            end
            if length(vtOpt)>1
                fprintf(model_fileID,'.param vt=%.12e\n',str2double(vt_opt));
            else
                fprintf(model_fileID,'.param vt=%.12e\n',vtOpt);                
            end
            fprintf(model_fileID,'.param etares=%.12e\n',str2double(etares_opt));
            fprintf(model_fileID,'.param vres=%.12e\n',str2double(vres_opt));
        end
        fprintf(model_fileID,'.param CH0=%.12e\n',CH0_opt);
        fprintf(model_fileID,'.param beta=%.12e\n',beta_opt);
        fprintf(model_fileID,'.param EI=0\n\n');
        fprintf(model_fileID,'.param gam=%.3e *gam=0 no SF\n', gam);
        fprintf(model_fileID,'.param gam0=0\n');
        fprintf(model_fileID,'.param isb=%.3e *isb=1 no SB\n', isb);
        fprintf(model_fileID,'.param RPP=1e9\n');
        fprintf(model_fileID,'.param ri=1\n\n');

        fprintf(model_fileID,'*Auxiliary functions\n');
        fprintf(model_fileID,'.param I0(x)=''V(H)>%.2e ? imin+(imax-imin)*limit(0,1,x) : 0''\n', H0_limit);
        fprintf(model_fileID,'.param A(x)=''amin+(alphamax-alphamin)*limit(0,1,x)''\n');
        fprintf(model_fileID,'.param RS(x)=''rsmin+(rsmax-rsmin)*limit(0,1,x)''\n');
        fprintf(model_fileID,'.param VSB(x)=''x>isb ? vt : vset''\n');
        fprintf(model_fileID,'.param ISF(x)=''gam==0 ? 1 : pow(limit(0,1,x),gam)-gam0''\n');
        fprintf(model_fileID,'.param S(x)=''1/(1+exp(-etaset*(x-VSB(I(RS)))))''\n');
        fprintf(model_fileID,'.param R(x)=''1/(1+exp(-etares*ISF(V(H))*(x-vres)))''\n');
        fprintf(model_fileID,'.param limit(y,z,x)=''x>z ? z : (x<y ? y : x)''\n\n');
        
        fprintf(model_fileID,'*Quasi-static model: H-V\n');
        fprintf(model_fileID,'GH   gnd! H 	cur=''min(R(V(C,n)),max(S(V(C,n)),V(H)))''\n');
        fprintf(model_fileID,'Rpar gnd! H 	R=1\n');
        fprintf(model_fileID,'CH   H    gnd! 	C=''CH0'' IC=''H0''\n\n');

        fprintf(model_fileID,'*Quasi-static model: I-V\n');
        fprintf(model_fileID,'RE   p    C 	R=''ri''\n');
        fprintf(model_fileID,'RS   C 	B 	R=''RS(V(H))''\n');
        fprintf(model_fileID,'GD   B    n 	cur=''I0(V(H))*(exp(beta*A(V(H))*V(B,n))-exp(-(1-beta)*A(V(H))*V(B,n)))''\n');
        fprintf(model_fileID,'RB   p    n 	R=''RPP''\n');
        fprintf(model_fileID,'GR   p    n 	cur=1e-12\n');     

        fprintf(model_fileID,'.ends memdiode\n');                 
    end
    fclose(model_fileID);    
end