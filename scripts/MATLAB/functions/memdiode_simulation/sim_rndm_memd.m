function [data_sim] = sim_rndm_memd(varargin)
    imaxRange=[1e-5,1e-2];
    iminRange=[1e-7,1e-5];
    alphaminRange=[1,5];
    alphamaxRange=[1,5];
    alphaRange=[1,5];
    rsRange=[1,100];
    rsminRange=[1,100];
    rsmaxRange=[1,50];
    etaresRange=[1,50];
    etasetRange=[1,100];
    vresRange=[-1.5 -0.1];
    vsetRange=[0.1 1.5];
    choRange=[1e-4 1e-1];
    beta=0.5;
    nRuns=10000;
    model='QMM';
    dbDistro='uniform';
    RUNLVL=5;
    freq=1;
    re_freq_vect=[5e2 2e2 1e3];
    Vapp_vect=[2];
    tstop=1/freq;
    tstep=tstop/200;
    tDELMAX=tstop/1000;
    dbFolder='/home/users/aguirref/nn_rs_uab/DBs/QMM';
    dbFile='QMM_sim_DB.mat';
    saveFlag=0;
    vt_mode='sin';
    useAlphaMaxMin=0;
    
    multiWaitbar( 'CloseAll' );
    
    if mod(nargin-0,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-0
            if strcmp(varargin{i},'imaxRange')
                imaxRange=varargin{i+1};
            end
            if strcmp(varargin{i},'iminRange')
                iminRange=varargin{i+1};
            end
            if strcmp(varargin{i},'etaresRange')
                etaresRange=varargin{i+1};
            end
            if strcmp(varargin{i},'vresRange')
                vresRange=varargin{i+1};
            end
            if strcmp(varargin{i},'etasetRange')
                etasetRange=varargin{i+1};
            end
            if strcmp(varargin{i},'vsetRange')
                vsetRange=varargin{i+1};
            end
            if strcmp(varargin{i},'betaOpt')
                betaRange=varargin{i+1};
            end
            if strcmp(varargin{i},'rsmaxRange')
                rsmaxRange=varargin{i+1};
            end
            if strcmp(varargin{i},'rsminRange')
                rsminRange=varargin{i+1};
            end
            if strcmp(varargin{i},'rsRange')
                rsRange=varargin{i+1};
            end
            if strcmp(varargin{i},'alphamaxRange')
                alphamaxRange=varargin{i+1};
                useAlphaMaxMin=1;
            end
            if strcmp(varargin{i},'alphaminRange')
                alphaminRange=varargin{i+1};
                useAlphaMaxMin=1;
            end
            if strcmp(varargin{i},'alphaRange')
                alphaRange=varargin{i+1};
            end
            if strcmp(varargin{i},'nRuns')
                nRuns=varargin{i+1};
            end
            if strcmp(varargin{i},'model')
                model=varargin{i+1};
            end
            if strcmp(varargin{i},'dbDistro')
                dbDistro=varargin{i+1};
            end
            if strcmp(varargin{i},'re_freq_vect')
                re_freq_vect=varargin{i+1};
            end
            if strcmp(varargin{i},'freq')
                freq=varargin{i+1};
            end
            if strcmp(varargin{i},'Vapp_vect')
                Vapp_vect=varargin{i+1};
            end
            if strcmp(varargin{i},'tstop')
                tstop=varargin{i+1};
            end
            if strcmp(varargin{i},'tstep')
                tstep=varargin{i+1};
            end
            if strcmp(varargin{i},'tDELMAX')
                tDELMAX=varargin{i+1};
            end     
            if strcmp(varargin{i},'RUNLVL')
                RUNLVL=varargin{i+1};
            end     
            if strcmp(varargin{i},'dbFile')
                dbFile=varargin{i+1};
            end     
            if strcmp(varargin{i},'vtMode')
                vt_mode=varargin{i+1};
            end                 
        end
    end   
    
    paramsRange.imax=imaxRange;
    paramsRange.imin=iminRange;
    if  useAlphaMaxMin==0
        paramsRange.alpha=alphaRange;
    else
        paramsRange.alphamin=alphaminRange;
        paramsRange.alphamax=alphamaxRange;
    end
    paramsRange.etaset=etasetRange;
    paramsRange.etares=etaresRange;
    paramsRange.vset=vsetRange;
    paramsRange.vres=vresRange;
    paramsRange.rs=rsRange;
    
    for Vapp_i=1:length(Vapp_vect)
        Vapp=Vapp_vect(Vapp_i);
        multiWaitbar('Vapp index...',Vapp_i/length(Vapp_vect));
        
        for re_freq_i=1:length(re_freq_vect)    
            re_freq=re_freq_vect(re_freq_i);
            multiWaitbar('re_freq index...',re_freq_i/length(re_freq_vect));
            
            if exist(fullfile(dbFolder,dbFile),'file')
        
                load(fullfile(dbFolder,dbFile));

                DB_QMM_index=(DB_QMM.Vapp_max==Vapp) &...
                             (DB_QMM.freq==freq) &...
                             (DB_QMM.re_freq==re_freq) &...
                             (DB_QMM.tstep==tstep) &...
                             (DB_QMM.tDELMAX==tDELMAX) &...
                             (any(string(DB_QMM.Model)==model,2)) &...
                             (any(string(DB_QMM.DB_distro)==dbDistro,2));
                         
                DB_QMM_index2search=find(DB_QMM_index);         
                
                DB_QMM_index2searchII=[];
                DB_QMM_index2searchII_counter=0;                             
                if ~isempty(DB_QMM_index2search)
                    for DB_QMM_index2search_i=1:length(DB_QMM_index2search)
                        if size(struct2table(DB_QMM.param_definitions{DB_QMM_index2search(DB_QMM_index2search_i),1}),2)==size(struct2table(paramsRange),2)
                            DB_QMM_index2searchII_counter=DB_QMM_index2searchII_counter+1;
                            DB_QMM_index2searchII(DB_QMM_index2searchII_counter)=DB_QMM_index2search(DB_QMM_index2search_i);
                        end
                    end
                end
                
                DB_QMM_index2searchIII=[];
                DB_QMM_index2searchIII_counter=0;
                if ~isempty(DB_QMM_index2searchII)
                    for DB_QMM_index2search_i=1:length(DB_QMM_index2searchII)
                        if ismember(struct2table(DB_QMM(find(DB_QMM_index),:).param_definitions{DB_QMM_index2searchII(DB_QMM_index2search_i),1}),struct2table(paramsRange))
                            DB_QMM_index2searchIII_counter=DB_QMM_index2searchIII_counter+1;
                            DB_QMM_index2searchIII(DB_QMM_index2searchIII_counter)=DB_QMM_index2searchII(DB_QMM_index2search_i);
                        end
                    end
                end
                
                if isempty(DB_QMM_index2searchIII)
                    index_DB_QMM=size(DB_QMM,1)+1;
                    start_nRUM=1;
                else
                    start_nRUM=size(DB_QMM.I_memd_array{DB_QMM_index2searchIII,1},1)+1; 
                    index_DB_QMM=DB_QMM_index2searchIII;
                end
            else
                start_nRUM=1;
                index_DB_QMM=1;
            end
                    
            for i=start_nRUM:nRuns
                
                saveFlag=1;

                multiWaitbar('Run number...',i/nRuns);
                
                if strcmpi(dbDistro,'uniform')
                    
                    imax=10^(log10(min(imaxRange))+(log10(max(imaxRange))-log10(min(imaxRange))).*rand(1,1));
                    imin=10^(log10(min(iminRange))+(log10(max(iminRange))-log10(min(iminRange))).*rand(1,1));
                    rs=10^(log10(min(rsRange))+(log10(max(rsRange))-log10(min(rsRange))).*rand(1,1));
                    etares=10^(log10(min(etaresRange))+(log10(max(etaresRange))-log10(min(etaresRange))).*rand(1,1));
                    etaset=10^(log10(min(etasetRange))+(log10(max(etasetRange))-log10(min(etasetRange))).*rand(1,1));
                    vres=min(vresRange)+(max(vresRange)-min(vresRange)).*rand(1,1);    
                    vset=min(vsetRange)+(max(vsetRange)-min(vsetRange)).*rand(1,1);    
                    alpha=min(alphaRange)+(max(alphaRange)-min(alphaRange)).*rand(1,1);    
                    alphamin=10^(log10(min(alphaminRange))+(log10(max(alphaminRange))-log10(min(alphaminRange))).*rand(1,1));    
                    alphamax=10^(log10(min(alphamaxRange))+(log10(max(alphamaxRange))-log10(min(alphamaxRange))).*rand(1,1));  
                    
                elseif strcmpi(dbDistro,'gaussian')
                    imin=1e3;
                    imax=abs(10^(normrnd(log10(imaxRange(1)),imaxRange(2))));
                    while imin>=imax
                        imin=abs(10^(normrnd(log10(iminRange(1)),iminRange(2))));
                    end
                    rs=abs(normrnd(rsRange(1),rsRange(2)));
                    etares=abs(normrnd(etaresRange(1),etaresRange(2)));
                    etaset=abs(normrnd(etasetRange(1),etasetRange(2)));
                    vres=-1*abs(normrnd(vresRange(1),vresRange(2)));
                    vset=abs(normrnd(vsetRange(1),vsetRange(2)));
                    if useAlphaMaxMin==0
                        alpha=abs(normrnd(alphaRange(1),alphaRange(2)));
                    else
                        alphamin=abs(normrnd(alphaminRange(1),alphaminRange(2)));
                        alphamax=abs(normrnd(alphamaxRange(1),alphamaxRange(2)));
                    end

                end

                if ~exist('DB_QMM','var')        
                    % The DB does not exist, and thereby it is created and the
                    % simulation results placed in the first entry of the table
                    if useAlphaMaxMin==0
                        sim_data=sim_memd_DB(beta, etaset, vset, etares, vres, imax, alpha, rs, imin, RUNLVL, tDELMAX, Vapp, tstep, tstop, re_freq, vt_mode);
                    else
                        sim_data=sim_memd_DB(beta, etaset, vset, etares, vres, imax, alpha, rs, imin, RUNLVL, tDELMAX, Vapp, tstep, tstop, re_freq, vt_mode, 'alphamin', alphamin, 'alphamax', alphamax);
                    end
                    
                    time=sim_data(:,1);
                    v_app=sim_data(:,2);
                    I_memd=sim_data(:,3);
                    v_h=sim_data(:,4);            
                    
                    if useAlphaMaxMin==0
                        QMM_params=table(imax,imin,alpha,etaset,etares,vset,vres,rs);
                        QMM_params.Properties.VariableNames = {'imax' 'imin' 'alpha' 'etaset' 'etares' 'vset' 'vres' 'rs'};    
                    else
                    	QMM_params=table(imax,imin,alphamax, alphamin,etaset,etares,vset,vres,rs);
                        QMM_params.Properties.VariableNames = {'imax' 'imin' 'alphamax' 'alphamin' 'etaset' 'etares' 'vset' 'vres' 'rs'};     
                    end

                    DB_QMM=table(Vapp, freq, re_freq, tstep, tDELMAX, {model}, {dbDistro}, {paramsRange}, {array2table(time')},{array2table(v_app')},{array2table(I_memd')},{array2table(v_h')},{QMM_params});
                    DB_QMM.Properties.VariableNames = {'Vapp_max' 'freq' 're_freq' 'tstep' 'tDELMAX' 'Model' 'DB_distro' 'param_definitions' 'time_array' 'Vapp_array' 'I_memd_array' 'V_h_array' 'QMM_params'};   
                    clear QMM_params

                else
                    % The gral. DB file exist, and an entry with a compatible set
                    % of general settings was found.            
                    if index_DB_QMM>size(DB_QMM,1)

                        if useAlphaMaxMin==0
                            sim_data=sim_memd_DB(beta, etaset, vset, etares, vres, imax, alpha, rs, imin, RUNLVL, tDELMAX, Vapp, tstep, tstop, re_freq, vt_mode);
                        else
                            sim_data=sim_memd_DB(beta, etaset, vset, etares, vres, imax, 3, rs, imin, RUNLVL, tDELMAX, Vapp, tstep, tstop, re_freq, vt_mode, 'alphamin', alphamin, 'alphamax', alphamax);
                        end
                        
                        time=sim_data(:,1);
                        v_app=sim_data(:,2);
                        I_memd=sim_data(:,3);
                        v_h=sim_data(:,4);            

                        if useAlphaMaxMin==0
                            QMM_params=table(imax,imin,alpha,etaset,etares,vset,vres,rs);
                            QMM_params.Properties.VariableNames = {'imax' 'imin' 'alpha' 'etaset' 'etares' 'vset' 'vres' 'rs'};    
                        else
                            QMM_params=table(imax,imin,alphamax, alphamin,etaset,etares,vset,vres,rs);
                            QMM_params.Properties.VariableNames = {'imax' 'imin' 'alphamax' 'alphamin' 'etaset' 'etares' 'vset' 'vres' 'rs'};     
                        end   

                        DB_QMM_aux=table(Vapp, freq, re_freq, tstep, tDELMAX, {model}, {dbDistro}, {paramsRange}, {array2table(time')},{array2table(v_app')},{array2table(I_memd')},{array2table(v_h')},{QMM_params});
                        DB_QMM_aux.Properties.VariableNames = {'Vapp_max' 'freq' 're_freq' 'tstep' 'tDELMAX' 'Model' 'DB_distro' 'param_definitions' 'time_array' 'Vapp_array' 'I_memd_array' 'V_h_array' 'QMM_params'};   
                        DB_QMM=[DB_QMM;DB_QMM_aux];
                        clear QMM_params DB_QMM_aux;     
                    else
                        if iscell(DB_QMM.QMM_params)
                            QMM_params_aux=DB_QMM.QMM_params{index_DB_QMM,1};
                        else
                            QMM_params_aux=DB_QMM.QMM_params;
                        end
                        
                        if useAlphaMaxMin==0
                            QMM_params_index=(QMM_params_aux.imax==imax) &...
                                               (QMM_params_aux.imin==imin) &...
                                               (QMM_params_aux.alpha==alpha) &...
                                               (QMM_params_aux.rs==rs) &...
                                               (QMM_params_aux.etares==etares) &...
                                               (QMM_params_aux.etaset==etaset) &...
                                               (QMM_params_aux.vset==vset) &...
                                               (QMM_params_aux.vres==vres);
                        else
                            QMM_params_index=(QMM_params_aux.imax==imax) &...
                                               (QMM_params_aux.imin==imin) &...
                                               (QMM_params_aux.alphamax==alphamax) &...
                                               (QMM_params_aux.alphamin==alphamin) &...
                                               (QMM_params_aux.rs==rs) &...
                                               (QMM_params_aux.etares==etares) &...
                                               (QMM_params_aux.etaset==etaset) &...
                                               (QMM_params_aux.vset==vset) &...
                                               (QMM_params_aux.vres==vres);                            
                        end
                        if any(QMM_params_index)==0

                            if useAlphaMaxMin==0
                                sim_data=sim_memd_DB(beta, etaset, vset, etares, vres, imax, alpha, rs, imin, RUNLVL, tDELMAX, Vapp, tstep, tstop, re_freq, vt_mode);
                            else
                                sim_data=sim_memd_DB(beta, etaset, vset, etares, vres, imax, 3, rs, imin, RUNLVL, tDELMAX, Vapp, tstep, tstop, re_freq, vt_mode, 'alphamin', alphamin, 'alphamax', alphamax);
                            end

                            time=sim_data(:,1);
                            v_app=sim_data(:,2);
                            I_memd=sim_data(:,3);
                            v_h=sim_data(:,4);   

                            DB_QMM.time_array{index_DB_QMM,1}=[DB_QMM.time_array{index_DB_QMM,1};array2table(time')];
                            DB_QMM.Vapp_array{index_DB_QMM,1}=[DB_QMM.Vapp_array{index_DB_QMM,1};array2table(v_app')];
                            DB_QMM.I_memd_array{index_DB_QMM,1}=[DB_QMM.I_memd_array{index_DB_QMM,1};array2table(I_memd')];
                            DB_QMM.V_h_array{index_DB_QMM,1}=[DB_QMM.V_h_array{index_DB_QMM,1};array2table(v_h')];

                            if useAlphaMaxMin==0
                                DB_QMM.QMM_params{index_DB_QMM,:}=[QMM_params_aux;table(imax,imin,alpha,etaset,etares,vset,vres,rs)];   
                            else
                                DB_QMM.QMM_params{index_DB_QMM,:}=[QMM_params_aux;table(imax,imin,alphamax,alphamin,etaset,etares,vset,vres,rs)];    
                            end                               
                            
                        end
                    end
                end
            end
            if saveFlag==1
                if ~exist(dbFolder,'dir')
                    mkdir(dbFolder);
                end
                save(fullfile(dbFolder,dbFile),'DB_QMM');
                saveFlag=0;
            end
        end
    end
    
    multiWaitbar( 'CloseAll' );
end