function [DATA] = procSim_MLP_results(dir_n_files,sim_settings,CPA_settings,sim_i,DATA_SIM_1,varargin)
    remap_opt='no';
    noise_opt='no';
    repeat_read=1;
    meas_mode='voltage';
    %% Optional values assignment
    if mod(nargin-5,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-5
            if strcmp(varargin{i},'remap')
                remap_opt=varargin{i+1};
            end
            if strcmp(varargin{i},'noise')
                noise_opt=varargin{i+1};
            end            
            if strcmp(varargin{i},'measMode')
                meas_mode=varargin{i+1};
            end            
        end
    end

    img_folder=dir_n_files.img_folder;
    global_figs=dir_n_files.global_figs;
    read_results_m=dir_n_files.read_results_m;
    read_results_m_remap=dir_n_files.read_results_m_remap;
    read_results_m_remap_fresh=dir_n_files.read_results_m_remap_fresh;
    read_results_m_noise=dir_n_files.read_results_m_noise;

    
    close_after_save=sim_settings.close_after_save;
    gen_plot=sim_settings.gen_plot;    
    time_offset=sim_settings.time_offset;
    input_vector_freq=sim_settings.input_vector_freq;
    in_vect_r_f_time=sim_settings.in_vect_r_f_time;
    connections=sim_settings.connections;
    correct_folder=sim_settings.correct_folder;
    labels_t10k=sim_settings.labels_t10k;
    map_str=sim_settings.map_str;
    %memdiode_model=sim_settings.memdiode_model;
    image_size=sim_settings.image_size;
    
    digits_order=CPA_settings.digits_order;
    partitions_N=CPA_settings.partitions_N;
    pos_neg=CPA_settings.pos_neg;
    series_resistance=CPA_settings.RS;
    R_cs=CPA_settings.R_cs;
    read_voltage=CPA_settings.read_voltage;
    write_voltage=CPA_settings.write_voltage;
    number_of_images=CPA_settings.number_of_images;
    neurons_per_layer=CPA_settings.neurons_per_layer;    
    neural_layers=CPA_settings.neural_layers;
    
    if strcmpi(noise_opt,'noisy')
        repeat_read=2;
    end
       
    for repeat_read_i=1:repeat_read
        if repeat_read_i==1
            if strcmpi(remap_opt,'faulty')
                list_files=dir(read_results_m_remap);
                for file_i=1:length(list_files)
                    if ~strcmpi(list_files(file_i,1).name,'..') && ~strcmpi(list_files(file_i,1).name,'.')
                        run(fullfile(list_files(file_i,1).folder,list_files(file_i,1).name));
                    end
                end                  
            elseif strcmpi(remap_opt,'fresh')
                list_files=dir(read_results_m_remap_fresh);
                for file_i=1:length(list_files)
                    if ~strcmpi(list_files(file_i,1).name,'..') && ~strcmpi(list_files(file_i,1).name,'.')
                        run(fullfile(list_files(file_i,1).folder,list_files(file_i,1).name));
                    end
                end
            else
                list_files=dir(read_results_m);
                for file_i=1:length(list_files)
                    if ~strcmpi(list_files(file_i,1).name,'..') && ~strcmpi(list_files(file_i,1).name,'.')
                        run(fullfile(list_files(file_i,1).folder,list_files(file_i,1).name));
                    end
                end
            end
        else
            list_files=dir(read_results_m_noise);
            for file_i=1:length(list_files)
                if ~strcmpi(list_files(file_i,1).name,'..') && ~strcmpi(list_files(file_i,1).name,'.')
                    run(fullfile(list_files(file_i,1).folder,list_files(file_i,1).name));
                end
            end            
        end
        
        V_neuron=[];
        for partitions_ii=1:partitions_N(end)
            for i=1:neurons_per_layer(end)/partitions_N(end)
                if exist(sprintf('v__neuron%s%d_p%d_1__', lower(neural_layers{end}), i, partitions_ii),'var')
                    str_rneuron=sprintf('V_neuron(%d,:)=v__neuron%s%d_p%d_1__;', i+(neurons_per_layer(end)/partitions_N(end))*(partitions_ii-1), lower(neural_layers{end}), i, partitions_ii);
                    eval(str_rneuron);
                end
            end
        end

        if isempty(V_neuron) || strcmpi(meas_mode,'current')
            for i=1:neurons_per_layer(end)
                str_rneuron=sprintf('Ineuron=');
                for pos_neg_i=1:2
                    str_rneuron=strcat(str_rneuron,'(');
                    for partitions_i=1:partitions_N(end-1)
                        for partitions_ii=1:partitions_N(end)
                            str_rneuron=strcat(str_rneuron,sprintf('i__rneuron%s%d_p%d_%d_%s__', lower(neural_layers{end}), i, partitions_i, partitions_ii, pos_neg{pos_neg_i}));
                            if partitions_i<partitions_N(end-1)
                                str_rneuron=strcat(str_rneuron,sprintf('+'));
                            end
                        end
                    end
                    str_rneuron=strcat(str_rneuron,')');
                    if pos_neg_i<2
                        str_rneuron=strcat(str_rneuron,sprintf('-'));
                    end
                end
                str_rneuron=strcat(str_rneuron,sprintf(';'));
                eval(str_rneuron);
                if i==1
                    eval(sprintf('I_neuron=Ineuron;'));
                else
                    eval(sprintf('I_neuron=vertcat(I_neuron,Ineuron);'));
                end
            end
        else
            I_neuron=V_neuron;
        end   
        % Acá entra solo si se realiza una simulación de ruido. primerpo
        % entra en laprimera opción para levantar la info de la señal sin
        % ruido, y despues en la segunda para procesar la infor de la señal
        % ruidosa.
        if repeat_read_i==1 && repeat_read==2
            I_neuron_noisless=I_neuron;
            TIME_noisless=TIME;
            signal_rms=rms(rms(I_neuron_noisless));
            clear I_neuron Ineuron clear ^i__rneuron ^v__neuron TIME
        elseif repeat_read_i==2 && repeat_read==2
            I_neuron_noisy=I_neuron;
            TIME_noisy=TIME;
            
            for neuron_i=1:neurons_per_layer(end) 
                [TIME_noisless_unique, index] = unique(TIME_noisless); 
                I_neuron_noisless_resampled(neuron_i,:)=interp1(TIME_noisless_unique,I_neuron_noisless(neuron_i,index),TIME_noisy);
            end
            
            noise_signals=I_neuron_noisy-I_neuron_noisless_resampled;
            noise_rms=rms(rms(noise_signals));
        end
    end
    
    if strcmpi(noise_opt,'noisy')
        I_neuron=avg_noisy_signal(TIME, I_neuron, sim_settings, 100);
    end

    [max_val,max_index]=max(I_neuron);

    I_neuron(neurons_per_layer(end)+1,:)=max_index;%-1;

    meas_time=time_offset+1/input_vector_freq/2+in_vect_r_f_time:1/input_vector_freq+in_vect_r_f_time:TIME(end);
    meas_time=vertcat(meas_time,zeros(4,size(meas_time,2)));

    num_classes=neurons_per_layer(end);
    for i=1:size(meas_time,2)
        aux_time=TIME-meas_time(1,i);
        min_index=min(find(aux_time>0));
        max_index=max(find(aux_time<0));
        if any(size(max_index)==0)
            max_index=length(aux_time);
        end
        if any(size(min_index)==0)
            min_index=length(aux_time);
        end        
        meas_time(2,i)=max_index;
        meas_time(3,i)=min_index;
        meas_time(4,i)=I_neuron(num_classes+1,min_index);
        meas_time(5,i)=I_neuron(num_classes+1,max_index);
    end

    %% Column re-mapping
    for i=1:size(meas_time,2)
        meas_time(4,i)=digits_order(meas_time(4,i),2);
        meas_time(4,i)=meas_time(4,i)-1;
        meas_time(5,i)=digits_order(meas_time(5,i),2);
        meas_time(5,i)=meas_time(5,i)-1;
    end
    %% label correct value:
    meas_time(6,:)=zeros(1,size(meas_time,2));
    meas_time(6,1:number_of_images)=labels_t10k(1:number_of_images);

    meas_time=meas_time(:,1:number_of_images);
    for i=0:(num_classes-1)
        DIGITS{i+1,1}=i;
        DIGITS{i+1,2}=find(meas_time(6,1:number_of_images)==i);

        ERRORS{i+1,1}=i;
        ERRORS{i+1,2}=find(meas_time(4,:)~=meas_time(6,:) & meas_time(6,:)==i);
        ERRORS{i+1,3}=sum(meas_time(4,:)~=meas_time(6,:) & meas_time(6,:)==i)/size(DIGITS{i+1,2},2);

        TP{i+1,1}=i;
        TP{i+1,2}=find(meas_time(4,:)==i & meas_time(6,:)==i);
        TP{i+1,3}=sum(meas_time(4,:)==i & meas_time(6,:)==i);
        
        TN{i+1,1}=i;
        TN{i+1,2}=find(meas_time(4,:)~=i & meas_time(6,:)~=i);
        TN{i+1,3}=sum(meas_time(4,:)~=i & meas_time(6,:)~=i);
        
        FP{i+1,1}=i;
        FP{i+1,2}=find(meas_time(4,:)==i & meas_time(6,:)~=i);
        FP{i+1,3}=sum(meas_time(4,:)==i & meas_time(6,:)~=i);
        
        FN{i+1,1}=i;
        FN{i+1,2}=find(meas_time(4,:)~=i & meas_time(6,:)==i);
        FN{i+1,3}=sum(meas_time(4,:)~=i & meas_time(6,:)==i);
        
        accuracy{i+1,1} = (TP{i+1,3}+TN{i+1,3})/(TP{i+1,3}+TN{i+1,3}+FP{i+1,3}+FN{i+1,3});
        precision{i+1,1} =TP{i+1,3} / (TP{i+1,3} + FP{i+1,3});
        sensitivity{i+1,1} = TP{i+1,3} / (TP{i+1,3} + FN{i+1,3}); 
        specificity{i+1,1} = TN{i+1,3} / (FP{i+1,3} + TN{i+1,3}); 
        F1_score{i+1,1} = 2*TP{i+1,3} /(2*TP{i+1,3} + FP{i+1,3} + FN{i+1,3});
        
        CURRENTS{i+1,1}=i;
        CURRENTS{i+1,2}=I_neuron(1:num_classes,meas_time(2,find(meas_time(6,:)==i)));
        CURRENTS{i+1,3}=mean(CURRENTS{i+1,2},2)';
        CURRENTS{i+1,4}=mean(CURRENTS{i+1,2}./sum(CURRENTS{i+1,2},1),2);
                               
    end

    FN_tot=0;
    FP_tot=0;
    TP_tot=0;
    TN_tot=0;
    for index_class=1:length(FP)  
        TP_tot=TP_tot+TP{index_class,3};
        TN_tot=TN_tot+TN{index_class,3};
        FP_tot=FP_tot+FP{index_class,3};
        FN_tot=FN_tot+FN{index_class,3};
    end

    PRECISION_gral = TP_tot/(FP_tot+TP_tot);
    SENSITIVITY_gral = TP_tot/(FN_tot+TP_tot);
    SPECIFICITY_gral = TN_tot/(FP_tot+TN_tot);
    F1_score_gral = 2*TP_tot/(2*TP_tot + FP_tot + FN_tot);
    ACCURACY_gral = (TP_tot+TN_tot)/(TP_tot+TN_tot+FP_tot+FN_tot);        

    comparison=vertcat(meas_time(6,:),meas_time(5,:));
    comparison(3,:)=comparison(1,:)-comparison(2,:);
    errors=find(comparison(3,:));
    comparison(4,:)=zeros(1,size(comparison,2));
    comparison(4,errors)=1;
    ERROR_gral=sum(comparison(4,:))/size(comparison,2);
    
    edges=-0.5:1:(num_classes-0.5);
    for digit=0:(num_classes-1)
        PROB(digit+1,:)=histcounts(meas_time(4,meas_time(6,:)==digit),edges)/size(DIGITS{digit+1,2},2);
        CURR_MAT(digit+1,:)=CURRENTS{digit+1,3};
        PROB_MAT(digit+1,:)=CURRENTS{digit+1,4};
    end

    if sim_i==1
        max_CURR=max(max(CURR_MAT))*1.05;
        min_CURR=min(min(CURR_MAT))*1.05;
    else
        min_CURR=DATA_SIM_1.min_max_CURR(1);
        max_CURR=DATA_SIM_1.min_max_CURR(2);
    end
        
    if gen_plot==1

        if ~exist(img_folder,'dir')
            mkdir(img_folder)
        end   

        if ~exist(global_figs,'dir')
            mkdir(global_figs);
        end

        barplots=figure('Name',sprintf('%s, %s, Rs=%.2f, Rcs=%.2f, N=%d, V_read=%f, V_write=%f, N_partitions=%d, %s', connections, correct_folder, series_resistance, R_cs, number_of_images, read_voltage, write_voltage, partitions_N, map_str),'Position',[0 0 1200 600]);

        subplot(1,2,1)
        h1=bar3(PROB);
        %h1=histogram2(PROB,'Normalization','probability');
        colorbar
        for k = 1:length(h1)
            zdata = h1(k).ZData;
            h1(k).CData = zdata;
            h1(k).FaceColor = 'interp';
        end
        zlabel('Digit probability per neuron');
        ylabel('Digit [#]');
        yticks([1:1:num_classes]);
        yticklabels(cellstr(split(num2str([0:1:num_classes-1]))));
        %yticklabels({'0','1','2','3','4','5','6','7','8','9'});
        xlabel({'Output','Neuron [#]'});
        xticks([1:1:num_classes]);
        xticklabels(cellstr(split(num2str([0:1:num_classes-1]))));
        %xticklabels({'0','1','2','3','4','5','6','7','8','9'});
        title({sprintf('%dx%d',image_size(1),image_size(2)), strrep(sprintf('%s',map_str),'_',' '), strrep(sprintf('%s',connections),'_',' '), strrep(sprintf('%s',correct_folder),'_',' '), sprintf('R_{s}=%.2g, R_{cs}=%.2g, N=%d, V_{read}=%.2g, V_{write}=%.2g, N_{part}=%d', series_resistance, R_cs, number_of_images, read_voltage, write_voltage, partitions_N)})

        axis([0 (num_classes+1) 0 (num_classes+1) 0 1.1]);

        subplot(1,2,2)
        h2=bar3(CURR_MAT);
        colorbar
        for k = 1:length(h2)
            zdata = h2(k).ZData;
            h2(k).CData = zdata;
            h2(k).FaceColor = 'interp';
        end
        zlabel({'Average current distribution','per neuron per digit [A]'})
        ylabel('Digit [#]');
        yticks([1:1:num_classes]);
        yticklabels(cellstr(split(num2str([0:1:num_classes-1]))));
        %yticklabels({'0','1','2','3','4','5','6','7','8','9'});
        xlabel({'Output','Neuron [#]'});
        xticks([1:1:num_classes]);
        xticklabels(cellstr(split(num2str([0:1:num_classes-1]))));
        %xticklabels({'0','1','2','3','4','5','6','7','8','9'});
        title({sprintf('%dx%d',image_size(1),image_size(2)), strrep(sprintf('%s',map_str),'_',' '), strrep(sprintf('%s',connections),'_',' '), strrep(sprintf('%s',correct_folder),'_',' '), sprintf('R_{s}=%.2g, R_{cs}=%.2g, N=%d, V_{read}=%.2g, V_{write}=%.2g, N_{part}=%d', series_resistance, R_cs, number_of_images, read_voltage, write_voltage, partitions_N)})

        if min_CURR == 0 && max_CURR == 0
            min_CURR=-1e-6;
            max_CURR=1e-6;
        end
        axis([0 11 0 11 min_CURR max_CURR]);
        caxis([min_CURR/1.25 max_CURR/1.25]);

        set(gcf,'color','white');

        print (gcf,'-dpng',fullfile(img_folder,strcat('stats.png')))
        savefig(fullfile(img_folder,strcat('stats.fig')));                         

        filename=fullfile(global_figs,'stats_vs_Rs.gif');

        frame=getframe(barplots);
        im=frame2im(frame);
        [imind,cm]=rgb2ind(im,256);

        if sim_i==1
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end

        confussionMatrix_plot=figure('Name',sprintf('%s, %s, Rs=%.2f, Rcs=%.2f, N=%d, V_read=%f, N_partitions=%d, %s', connections, correct_folder, series_resistance, R_cs, number_of_images, read_voltage, partitions_N, map_str),'Position',[0 0 1200 600]);


        subplot(1,2,1)
        h1_=imagesc(PROB);
        colorbar

        zlabel('Digit probability per neuron');
        ylabel('Digit [#]');
        yticks([1:1:num_classes]);
        yticklabels(cellstr(split(num2str([0:1:num_classes-1]))));
        %yticklabels({'0','1','2','3','4','5','6','7','8','9'});
        xlabel({'Output','Neuron [#]'});
        xticks([1:1:num_classes]);
        xticklabels(cellstr(split(num2str([0:1:num_classes-1]))));
        %xticklabels({'0','1','2','3','4','5','6','7','8','9'});
        title({sprintf('%dx%d',image_size(1),image_size(2)), strrep(sprintf('%s',map_str),'_',' '), strrep(sprintf('%s',connections),'_',' '), strrep(sprintf('%s',correct_folder),'_',' '), sprintf('R_{s}=%.2g, R_{cs}=%.2g, N=%d, V_{read}=%.2g, N_{part}=%d', series_resistance, R_cs, number_of_images, read_voltage, partitions_N)})
        axis square
        axis([0 num_classes+1 0 num_classes+1]);

        subplot(1,2,2)
        h2_=imagesc(CURR_MAT);
        colorbar

        zlabel({'Average current distribution','per neuron per digit [A]'})
        ylabel('Digit [#]');
        yticks([1:1:num_classes]);
        yticklabels(cellstr(split(num2str([0:1:num_classes-1]))));
        %yticklabels({'0','1','2','3','4','5','6','7','8','9'});
        xlabel({'Output','Neuron [#]'});
        xticks([1:1:num_classes]);
        xticklabels(cellstr(split(num2str([0:1:num_classes-1]))));
        %xticklabels({'0','1','2','3','4','5','6','7','8','9'});
        title({sprintf('%dx%d',image_size(1),image_size(2)), strrep(sprintf('%s',map_str),'_',' '), strrep(sprintf('%s',connections),'_',' '), strrep(sprintf('%s',correct_folder),'_',' '), sprintf('R_{s}=%.2g, R_{cs}=%.2g, N=%d, V_{read}=%.2g, N_{part}=%d', series_resistance, R_cs, number_of_images, read_voltage, partitions_N)})

        axis square
        axis([0 num_classes+1 0 num_classes+1]);
        caxis([min_CURR/1.25 max_CURR/1.25]);

        set(gcf,'color','white');

        print (gcf,'-dpng',fullfile(img_folder,strcat('stats_map.png')))
        savefig(fullfile(img_folder,strcat('stats_map.fig'))); 

        filename=fullfile(global_figs,'stats_vs_Rs_map.gif');

        frame=getframe(confussionMatrix_plot);
        im=frame2im(frame);
        [imind,cm]=rgb2ind(im,256);

        if sim_i==1
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end

    DATA.Vr=read_voltage;
    DATA.Vw=write_voltage;
    DATA.Rs=series_resistance;
    DATA.R_cs=R_cs;
    DATA.N=number_of_images;
    DATA.CURRENTS=CURRENTS;
    DATA.ERRORS=ERRORS;
    DATA.FP=FP;
    DATA.FN=FN;
    DATA.TP=TP;
    DATA.TN=TN;
    DATA.ACCURACY=accuracy;
    DATA.PRECISION=precision;
    DATA.SENSITIVITY=sensitivity; 
    DATA.SPECIFICITY=specificity; 
    DATA.F1_SCORE=F1_score;
    DATA.PRECISION_gral = PRECISION_gral;
    DATA.SENSITIVITY_gral = SENSITIVITY_gral;
    DATA.SPECIFICITY_gral = SPECIFICITY_gral;
    DATA.F1_score_gral = F1_score_gral;
    DATA.ACCURACY_gral = ACCURACY_gral;        
    DATA.DIGITS=DIGITS;
    DATA.CURR_MAT=CURR_MAT;
    DATA.PROB=PROB;
    DATA.PROB_MAT=PROB_MAT;
    DATA.ERROR_gral=ERROR_gral;
    DATA.imgSize=image_size;
    DATA.min_max_CURR=[min_CURR max_CURR];
    if strcmpi(noise_opt,'noisy')
    DATA.SIGNAL_RMS=signal_rms;
    DATA.NOISE_RMS=noise_rms;
    DATA.SNR=signal_rms/noise_rms;
    end
    
    if close_after_save==1
        close(confussionMatrix_plot);
        close(barplots);
    end
end

