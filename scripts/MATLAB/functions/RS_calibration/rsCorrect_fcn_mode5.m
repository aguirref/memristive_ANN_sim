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
    %directory2='/home/Proyectos/Proyecto_NEUROMORPHICS_UTN-UAB/docs';
    
    calibrate='YES';
    
    crit_max=1;
    crit_min=1e-4;
    
    accuracy_max=0;
    
    for np_i=1:length(partitions)-1
        partitions_aux(np_i,1)=partitions(np_i);
        partitions_aux(np_i,2)=partitions(np_i+1);
    end
    
    partitions_aux2=partitions;
    partitions=partitions_aux;
    
    calibrate_loop='yes';
    
    loop_cnt=1;
    level=1;
    
    multiWaitbar( 'CloseAll' );
        
    calc_results_mat=[0 0];
    
    nPoints=5;
    fig_aux=figure();
    while strcmpi(calibrate_loop,'yes')
        
        %multiWaitbar('Calibrating CPA',(cal_criterion_i-1)/(length(cal_criterion_vec)));
        cal_criterion_vec=logspace(log10(crit_min),log10(crit_max),5);
        
        for cal_criterion_i=1:length(cal_criterion_vec)

            cal_criterion=cal_criterion_vec(cal_criterion_i);
            
            if ~any(ismember(calc_results_mat(:,1),cal_criterion))
            
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
                    %
                end

                load(fullfile(directory2,'G_cal.mat'));

                if loop_cnt==1
                    cal_results{cal_criterion_i,1}=cal_criterion;
                    cal_results{cal_criterion_i,2}=G_norm_aux;
                    cal_results{cal_criterion_i,3}=Cal_accuracy;
                else
                    inserted_row=0;
                    for cal_results_j=1:size(cal_results,1)-1
                        if cal_criterion>cal_results{cal_results_j,1} && cal_criterion<cal_results{cal_results_j+1,1}
                            cal_results(cal_results_j+1:end+1,:)=cal_results(cal_results_j:end,:);

                            cal_results{cal_results_j+1,1}=cal_criterion;
                            cal_results{cal_results_j+1,2}=G_norm_aux;
                            cal_results{cal_results_j+1,3}=Cal_accuracy;
                            
                            inserted_row=1;
                            
                            break
                        end
                    end
                    if inserted_row==0 && cal_criterion>cal_results{cal_results_j,1}
                        cal_results{cal_results_j+2,1}=cal_criterion;
                        cal_results{cal_results_j+2,2}=G_norm_aux;
                        cal_results{cal_results_j+2,3}=Cal_accuracy;
                        
                    elseif inserted_row==0 && cal_criterion<cal_results{1,1}
                        
                        cal_results(2:end+1,:)=cal_results(1:end,:);
                        
                        cal_results{1,1}=cal_criterion;
                        cal_results{1,2}=G_norm_aux;
                        cal_results{1,3}=Cal_accuracy;
                    end
                end
            end
        
        end
        
        calc_results_mat=cell2mat(cal_results(:,[1 3]));
        semilogx(calc_results_mat(:,1),calc_results_mat(:,2),'Marker','none')
        
        loop_cnt=loop_cnt+1;
        
        [~,max_location]=max(calc_results_mat(:,2));
        maximum_points=find(calc_results_mat(:,2)==calc_results_mat(max_location,2));
        if length(maximum_points)>1
            max_location=maximum_points(end);
        end

        [~,idxsIntoA]=intersect(calc_results_mat(:,1),cal_criterion_vec','stable');
        if length(unique(calc_results_mat(idxsIntoA,2)))==1 || (crit_max-crit_min)<1e-4
            calibrate_loop='no';    
            G_norm=cal_results{max(idxsIntoA),2};
            final_calibration_factor=cal_results{max(idxsIntoA),1};
        elseif max_location>1 && max_location<size(calc_results_mat,1)
            crit_max=calc_results_mat(max_location+1);
            crit_min=calc_results_mat(max_location-1);
            level=level+1;
            
            cal_criterion_vec=logspace(log10(crit_min),log10(crit_max),5);
        elseif length(maximum_points)==1 && max_location==1 &&  level==1
            nPoints=nPoints+1;
            crit_min=crit_min/10;
            cal_criterion_vec=logspace(log10(crit_min),log10(crit_max),nPoints);
        elseif length(maximum_points)==1 && max_location==size(calc_results_mat,1) && level==1
            nPoints=nPoints+1;
            crit_max=crit_max*10;
            cal_criterion_vec=logspace(log10(crit_min),log10(crit_max),nPoints);
        elseif length(maximum_points)>1 
            max_location=maximum_points(1);
            crit_max=calc_results_mat(max_location+1);
            crit_min=calc_results_mat(max_location-1);
            level=level+1;
        end
        
        
        
        %multiWaitbar('Calibrating CPA',(cal_criterion_i)/(length(cal_criterion_vec)))
    end
         
    for j=1:length(G_norm)
        G_norm{j,1}=G_norm{j,1}';
    end
    
    partitions=partitions_aux2;
    
    %multiWaitbar( 'CloseAll' );
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

    save(fullfile(rs_corrected_weights,'lambda_vals.mat'),'W_init','W_matrix','W_init_pos','W_init_neg','W_matrix_pos','W_matrix_neg','cal_results','final_calibration_factor');    
    
    if strcmpi(mode,'log1') || strcmpi(mode,'log2') || strcmpi(mode,'log3') || strcmpi(mode,'log4')

    elseif strcmpi(mode,'lin1') || strcmpi(mode,'lin2') || strcmpi(mode,'lin3') || strcmpi(mode,'lin4')

    elseif strcmpi(mode,'memdiode1') || strcmpi(mode,'memdiode2') || strcmpi(mode,'memdiode3') || strcmpi(mode,'memdiode4')
       
    end
    
    corrected_WM.W_init_pos=W_init_pos;
    corrected_WM.W_init_neg=W_init_neg;
    close(fig_aux);
end

