function corrected_WM = rsCorrect_fcn_mode3(CPA_settings, sim_settings, dir_n_files, G, mode,digits_order,I_min, I_max,model_ver,memdiode_model,V_read,Rs,partitions)
%RSCORRECT_FCN_MODE2 Summary of this function goes here
%   Detailed explanation goes here

    changes=[];
    prj_dir=fullfile('/home/users','aguirref','nn_rs_uab');
    
    images_t10k=CPA_settings.images_t10k;
    labels_t10k=CPA_settings.labels_t10k;
    neural_layers=CPA_settings.neural_layers;
    H_min=sim_settings.map_array{1,1};
    H_max=sim_settings.map_array{1,2};
    R_shunt=CPA_settings.R_shunt;
    rs_corrected_weights=dir_n_files.rs_corrected_weights;
    number_of_images=sim_settings.number_of_images;
    init_cal_criterion=CPA_settings.init_cal_criterion;
    
    for i=1:size(digits_order,1)
        if size(unique(digits_order(i,:)),2)>1
            if ~isempty(changes)
                for ii=1:size(changes,1)
                    if unique(changes(ii,:))~=unique(digits_order(i,:))
                       changes=vertcat(changes,digits_order(i,:));
                    end
                end
            else
                changes=vertcat(changes,digits_order(i,:));
            end
        end
    end
   
    if ~isempty(changes)
        for i=size(changes,1)
            sour_col=G(:,changes(i,1));
            dest_col=G(:,changes(i,2));

            G(:,changes(i,2))=sour_col;
            G(:,changes(i,1))=dest_col;
        end
    end
    
    if ~exist(rs_corrected_weights,'dir')
        mkdir(rs_corrected_weights);
    end
    
    RHRS=V_read/I_min;
    RLRS=V_read/I_max;
       
    output_side='up';
    
    accuracy_calc='YES';
            
    cant_muestras=number_of_images;
    
    for j=1:length(G)
        G_trasp{j,1}=G{j,1}';
    end
    
    %directory='/home/users/aguirref/nn_rs_uab/scripts/Python/calibration_tool';
    directory='/home/Proyectos/Proyecto_NEUROMORPHICS_UTN-UAB/docs';
    directory2=rs_corrected_weights;
    
    calibrate='YES';
    
    cal_criterion_vec=logspace(log10(5e-3),log10(5e-1),15);
    
    accuracy_max=0;
    
    for np_i=1:length(partitions)-1
        partitions_aux(np_i,1)=partitions(np_i);
        partitions_aux(np_i,2)=partitions(np_i+1);
    end
    
    partitions_aux2=partitions;
    partitions=partitions_aux;
    for cal_criterion_i=1:length(cal_criterion_vec)
        
        cal_criterion=cal_criterion_vec(cal_criterion_i);
        
        delete(fullfile(directory2,'G_cal.mat'));
        delete('G_cal.mat');
        delete(fullfile(directory2,'G_pos.csv'));
        delete('G_pos.csv');
        delete(fullfile(directory2,'G_neg.csv'));
        delete('G_neg.csv');
        delete(fullfile(directory2,'compensation_variables.mat'));
        delete('compensation_variables.mat');

        save(fullfile(directory2,'compensation_variables.mat'),'G_trasp','RHRS','RLRS','Rs','accuracy_calc','calibrate','cal_criterion','V_read','cant_muestras','partitions','images_t10k','labels_t10k','output_side','directory');

        system_command = sprintf('python36 %s/calibration.py ''''%s''''',directory,directory2);
        [status,cmdout]=system(system_command,'-echo');

        
        for j=1:length(G)
            loaded_G_cal_pos=load(fullfile(directory2,sprintf('G_pos%d.csv',j)));
            loaded_G_cal_neg=load(fullfile(directory2,sprintf('G_neg%d.csv',j)));

            G_norm_aux{j,1}=loaded_G_cal_pos-loaded_G_cal_neg;

            load(fullfile(directory2,'G_cal.mat'));
        end
        
        if Cal_accuracy>=accuracy_max
            accuracy_max=Cal_accuracy;
            criterion_index_max=cal_criterion_i;
            G_norm=G_norm_aux;
        end
        
        cal_results{cal_criterion_i,1}=cal_criterion;
        cal_results{cal_criterion_i,2}=G_norm_aux;
        cal_results{cal_criterion_i,3}=Cal_accuracy;

    end
         
    for j=1:length(G_norm)
        G_norm{j,1}=G_norm{j,1}';
    end
    
    partitions=partitions_aux2;
    
    multiWaitbar( 'CloseAll' );
    for layer_i=1:length(neural_layers)-1
        multiWaitbar('Calculating CPA weights',(layer_i-1)/(length(neural_layers)-1));

        G_norm_aux=G_norm{layer_i,1};

        [W_init_aux,W_matrix_aux]=G2W_map(G_norm_aux,H_max,H_min,mode,digits_order, I_min, I_max, model_ver,memdiode_model,V_read, 0, 0, 'layer_i',layer_i,'total_layers',length(neural_layers)-1);

        W_init_pos_aux=W_init_aux(:,:,1);
        W_init_neg_aux=W_init_aux(:,:,2);

        W_matrix_pos_aux=W_matrix_aux(:,:,1)*R_shunt;
        W_matrix_neg_aux=W_matrix_aux(:,:,2)*R_shunt;

        W_init{layer_i,1}=W_init_aux;
        W_matrix{layer_i,1}=W_matrix_aux;

        W_init_pos{layer_i,1}=W_init_pos_aux;
        W_init_neg{layer_i,1}=W_init_neg_aux;
        W_matrix_pos{layer_i,1}=W_matrix_pos_aux;
        W_matrix_neg{layer_i,1}=W_matrix_neg_aux;

        multiWaitbar('Calculating CPA weights',layer_i/(length(neural_layers)-1));
    end

    save(fullfile(rs_corrected_weights,'lambda_vals.mat'),'W_init','W_matrix','W_init_pos','W_init_neg','W_matrix_pos','W_matrix_neg','cal_results')
    
    
    if strcmpi(mode,'log1') || strcmpi(mode,'log2') || strcmpi(mode,'log3') || strcmpi(mode,'log4')

    elseif strcmpi(mode,'lin1') || strcmpi(mode,'lin2') || strcmpi(mode,'lin3') || strcmpi(mode,'lin4')

    elseif strcmpi(mode,'memdiode1') || strcmpi(mode,'memdiode2') || strcmpi(mode,'memdiode3') || strcmpi(mode,'memdiode4')
       
    end
    
    corrected_WM.W_init_pos=W_init_pos;
    corrected_WM.W_init_neg=W_init_neg;

end

