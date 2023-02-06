function [NN_elements]=train_MLP_fcn(resol, limit_output, in_polarity, database, neu_per_layer, varargin)
    valChecks=15;
    save_pref=1;
    all_logsig=0;
    trainFcn = 'trainscg';                          % use scaled conjugate gradient for training

    if mod(nargin-5,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-5
            if strcmp(varargin{i},'save?')
                save_pref=varargin{i+1};
            end
            if strcmp(varargin{i},'valChecks')
                valChecks=varargin{i+1};
            end
            if strcmp(varargin{i},'all_logsig')
                all_logsig=varargin{i+1};
            end
            if strcmp(varargin{i},'trainFcn')
                trainFcn=varargin{i+1};
            end
            if strcmp(varargin{i},'mcFolder')
               train_MC_folder=varargin{i+1};
            end            
            if strcmp(varargin{i},'train_tool')
               train_tool=varargin{i+1};
            end   
            if strcmp(varargin{i},'ANN_type')
               ANN_type=varargin{i+1};
            end   
             if strcmp(varargin{i},'learning_rate')
               learning_rate=varargin{i+1};
            end   
            if strcmp(varargin{i},'num_Epochs')
               num_Epochs=varargin{i+1};
            end 
            if strcmp(varargin{i},'quantization')
               quantization=varargin{i+1};
            end                
        end
    end
    
    if limit_output==1 && strcmpi(in_polarity,'unipolar')
        lims='limited_range';
        polarity=0;
    elseif limit_output==0  && strcmpi(in_polarity,'unipolar')
        lims='unlimited_range';
        polarity=0;
    elseif limit_output==1  && strcmpi(in_polarity,'bipolar')
        lims='unlimited_range_bipolar';
        polarity=0.5;
    elseif limit_output==0  && strcmpi(in_polarity,'bipolar')
        lims='unlimited_range_bipolar';
        polarity=0.5;
    end
    
    if strcmpi(trainFcn,'trainscg')
        learning_algorithm='Scaled_Conjugate_Gradient';
    else
        learning_algorithm='Scaled_Conjugate_Gradient';
    end
    
    if all_logsig==1
        learning_algorithm=strcat(learning_algorithm,'_allLogsig');
    end

    mat_file_DB=fullfile('..','..','..','DBs',database,lims,sprintf('%dby%d',resol(1), resol(2)),sprintf('train-images_%dx%d.mat',resol(1),resol(2)));
    if strcmpi(train_tool,'MATLAB')
        mat_file_weights = fullfile('..','..','..','DBs',database, lims, sprintf('%dby%d',resol(1), resol(2)), sprintf('%d_layers',length(neu_per_layer)), strjoin(string(neu_per_layer),'_'), train_tool, ANN_type, learning_algorithm, train_MC_folder);
    else
        if strcmpi(ANN_type,'DNN')
            if ismissing(quantization)
                mat_file_weights = fullfile('..','..','..','DBs',database, lims, sprintf('%dby%d',resol(1), resol(2)), sprintf('%d_layers',length(neu_per_layer)), strjoin(string(neu_per_layer),'_'), train_tool, ANN_type, learning_algorithm, train_MC_folder);
            else
                mat_file_weights = fullfile('..','..','..','DBs',database, lims, sprintf('%dby%d',resol(1), resol(2)), sprintf('%d_layers',length(neu_per_layer)), strjoin(string(neu_per_layer),'_'), train_tool, ANN_type, sprintf('%d_levels', quantization), learning_algorithm, train_MC_folder);                                        
            end
        else
            mat_file_weights = fullfile('..','..','..','DBs',database, lims, sprintf('%dby%d',resol(1), resol(2)), sprintf('%d_layers',length(neu_per_layer)), strjoin(string(neu_per_layer),'_'), train_tool, ANN_type, learning_algorithm, train_MC_folder);
        end
    end

    if exist(mat_file_DB,'file')
        load(mat_file_DB)
    else
        if strcmpi(database,'MNIST') || strcmpi(database,'F-MNIST') || strcmpi(database,'K-MNIST')
            ref_resol=[28 28];

            % and this one unzip them
            if ~exist(fullfile('..','..','..','DBs',database,'train-images-idx3-ubyte'),'file')
                gunzip(fullfile('..','..','..','DBs',database,'train-images-idx3-ubyte.gz'));
            end
            if ~exist(fullfile('..','..','..','DBs',database,'train-labels-idx1-ubyte'),'file')
                gunzip(fullfile('..','..','..','DBs',database,'train-labels-idx1-ubyte.gz'));
            end
            if ~exist(fullfile('..','..','..','DBs',database,'t10k-images-idx3-ubyte'),'file')
                gunzip(fullfile('..','..','..','DBs',database,'t10k-images-idx3-ubyte.gz'));
            end
            if ~exist(fullfile('..','..','..','DBs',database,'t10k-labels-idx1-ubyte'),'file')
                gunzip(fullfile('..','..','..','DBs',database,'t10k-labels-idx1-ubyte.gz'));
            end 

            images = loadMNISTImages(fullfile('..','..','..','DBs',database,'train-images-idx3-ubyte'));
            labels = loadMNISTLabels(fullfile('..','..','..','DBs',database,'train-labels-idx1-ubyte'))';

            images_t10k = loadMNISTImages(fullfile('..','..','..','DBs',database,'t10k-images-idx3-ubyte'));
            labels_t10k = loadMNISTLabels(fullfile('..','..','..','DBs',database,'t10k-labels-idx1-ubyte'))';        
        elseif strcmpi(database,'CIFAR-10') || strcmpi(database,'SVHN') || strcmpi(database,'YALE-FACE') || strcmpi(database,'YALE-FACE-B') || strcmpi(database,'ORL')
            ref_resol=[32 32];

            images = load(fullfile('..','..','..','DBs',database,'train-images-idx3-ubyte'),'data_grayscale');
            images = images.data_grayscale;
            labels = load(fullfile('..','..','..','DBs',database,'train-labels-idx1-ubyte'),'data_labels')';
            labels = double(labels.data_labels);

            images_t10k = load(fullfile('..','..','..','DBs',database,'t10k-images-idx3-ubyte'),'data_test_grayscale');
            images_t10k = images_t10k.data_test_grayscale;
            labels_t10k = load(fullfile('..','..','..','DBs',database,'t10k-labels-idx1-ubyte'),'data_test_labels')';  
            labels_t10k = double(labels_t10k.data_test_labels);

        end

        for i = 1:size(images,2)                                        % preview first 36 samples
            digit = reshape(images(:, i), [ref_resol(1), ref_resol(2)]);                     % row = 28 x 28 image
            digit_resized = imresize(digit,[resol(1) resol(2)]);
            if limit_output==1
                digit_resized(digit_resized>1)=1;
                digit_resized(digit_resized<0)=0;
            end
            eval(sprintf('images_%dx%d(:,i) = reshape(digit_resized, [resol(1)*resol(2),1]);',resol(1),resol(2)));
        end
        eval(sprintf('labels_%dx%d = labels;',resol(1),resol(2)));

        for i = 1:size(images_t10k,2)                                        % preview first 36 samples
            digit = reshape(images_t10k(:, i), [ref_resol(1), ref_resol(2)]);                     % row = 28 x 28 image
            digit_resized = imresize(digit,[resol(1) resol(2)]);
            if limit_output==1
                digit_resized(digit_resized>1)=1;
                digit_resized(digit_resized<0)=0;
            end
            eval(sprintf('images_t10k_%dx%d(:,i) = reshape(digit_resized, [resol(1)*resol(2),1]);',resol(1),resol(2)));
        end
        eval(sprintf('labels_t10k_%dx%d = labels_t10k;',resol(1),resol(2)));

        if ~exist(fullfile('..','..','..','DBs',database, lims, sprintf('%dby%d',resol(1), resol(2))),'dir')
            mkdir(fullfile('..','..','..','DBs',database, lims, sprintf('%dby%d',resol(1), resol(2))));
        end

        if save_pref == 1
            save(mat_file_DB,...
                                sprintf('images_%dx%d',resol(1),resol(2)),...
                                sprintf('labels_%dx%d',resol(1),resol(2)),...
                                sprintf('images_t10k_%dx%d',resol(1),resol(2)),...
                                sprintf('labels_t10k_%dx%d',resol(1),resol(2)));
        end
    end
    eval(sprintf('x = images_%dx%d-polarity;',resol(1),resol(2)));
    eval(sprintf('t = labels_%dx%d;',resol(1),resol(2)));

    if strcmpi(train_tool,'Python')

        system_command=sprintf(['''import sys\n' ...
                                'sys.path.append("%s")\n' ...
                                'import train_ANN\n' ...
                                'train_ANN.train_ANN("%s", "%s", [%s], "%s", 0, %d, %s, %d)'''],fullfile('..','..','Python','train_ANN'),mat_file_DB, mat_file_weights, strrep(strrep(strrep(mat2str(neu_per_layer(2:end)),' ',','),']',''),'[',''), ANN_type, num_Epochs, learning_rate, quantization);
        system(sprintf('python3 -c %s', system_command));
        load(fullfile(mat_file_weights,'G_values.mat'),'G_real','G_pos','G_neg');
        if ~iscell(G_neg)
            G_neg={squeeze(G_neg)'};
        elseif size(G_neg,1)==1
            G_neg=G_neg';
            G_neg = cellfun(@transpose,G_neg,'UniformOutput',false);
        end
        if ~iscell(G_pos)
            G_pos={squeeze(G_pos)'};
        elseif size(G_pos,1)==1
            G_pos=G_pos';
            G_pos = cellfun(@transpose,G_pos,'UniformOutput',false);
        end
        if ~iscell(G_real)
            G_real={squeeze(G_real)'};
        elseif size(G_real,1)==1
            G_real=G_real';
            G_real = cellfun(@transpose,G_real,'UniformOutput',false);
        end
        save(fullfile(mat_file_weights,'G_values.mat'),'G_real','G_neg','G_pos');
    else
        A=patternnet(neu_per_layer(2:end),trainFcn);
        A.inputs{1,1}.size=resol(1)*resol(2);
        A.layers{length(neu_per_layer),1}.size=neu_per_layer(end);
        A.LW{length(neu_per_layer),length(neu_per_layer)-1}=diag(ones(neu_per_layer(end),1));
        A.layerweights{length(neu_per_layer),length(neu_per_layer)-1}.learn=0;
        A.inputweights{1,1}.learnFcn='learncon';
        for layer_counter_i=1:length(neu_per_layer)
            if layer_counter_i<length(neu_per_layer)
                if layer_counter_i<length(neu_per_layer)-1
                    if strcmpi(in_polarity,'unipolar')
                        A.layers{layer_counter_i,1}.transferFcn='logsig';  
                    elseif strcmpi(in_polarity,'bipolar')
                        A.layers{layer_counter_i,1}.transferFcn='tansig'; 
                    else
                        A.layers{layer_counter_i,1}.transferFcn='logsig'; 
                    end
                else
                    if all_logsig==0
                        A.layers{layer_counter_i,1}.transferFcn='purelin';  
                    else
                        if strcmpi(in_polarity,'unipolar')
                            A.layers{layer_counter_i,1}.transferFcn='logsig';  
                        elseif strcmpi(in_polarity,'bipolar')
                            A.layers{layer_counter_i,1}.transferFcn='tansig'; 
                        else
                            A.layers{layer_counter_i,1}.transferFcn='logsig'; 
                        end
                    end
                end
            end
            A.biases{layer_counter_i,1}.learn=0;
            A.biases{layer_counter_i,1}.learn=0;
        end
        A.biasConnect=zeros(length(neu_per_layer),1);
    
        A.divideFcn='dividerand';
        A.divideMode='sample';
        A.divideParam.trainRatio=80/100;
        A.divideParam.valRatio=20/100;
        A.divideParam.testRatio=0/100;
    
        A.trainParam.max_fail=valChecks;
        A.trainParam.min_grad=1e-7;
        A.trainParam.epochs=5000;
    
        num_classes=neu_per_layer(end);
        
        t_vec=full(ind2vec(t+1,num_classes));
    
        %profile clear
        %profile on -timer 'cpu' 
        tStart = cputime;
        [A,tr_a] = train(A,x,t_vec);
        tStop = cputime - tStart;
        %p=profile('info');
        for layer_i=1:length(neu_per_layer)-1
            if layer_i==1
                G_real_aux=A.IW{1,1};
            else
                G_real_aux=A.LW{layer_i,layer_i-1};
            end
            
            G_aux=G_real_aux;
            G_aux(G_aux<0)=0;
    
            G_enf_aux=G_aux;
            G_enf_aux=log(G_enf_aux+1)/log(2);
    
            G_pos_aux=G_real_aux;
            G_pos_aux(G_pos_aux<0)=0;
    
            G_neg_aux=G_real_aux;
            G_neg_aux(G_neg_aux>0)=0;
    
            G_neg_aux=G_neg_aux*-1;
    
            G_bin_aux=G_pos_aux;
            G_bin_aux(G_bin_aux>0)=1;
    
            G_bin2_aux=G_pos_aux;
            G_bin2_aux(G_bin2_aux>0.5)=1;
            G_bin2_aux(G_bin2_aux<0.5)=0;   
            
            G_real{layer_i,1}=G_real_aux;
            G{layer_i,1}=G_aux;
            G_enf{layer_i,1}=G_enf_aux;
            G_pos{layer_i,1}=G_pos_aux;
            G_neg{layer_i,1}=G_neg_aux;
            G_bin{layer_i,1}=G_bin_aux;
            G_bin2{layer_i,1}=G_bin2_aux;        
        end
     
        eval(sprintf('results=A(images_t10k_%dx%d-polarity);',resol(1),resol(2)));
        [~,results]=max(results);
        results=results-1;
        eval(sprintf('comparison=(results==labels_t10k_%dx%d);',resol(1),resol(2)));
        Accuracy=sum(comparison)/length(comparison)*100;
        
        eval(sprintf('images_test = (images_t10k_%dx%d-polarity);',resol(1),resol(2)));
        eval(sprintf('labels_test = labels_t10k_%dx%d;',resol(1),resol(2)));
        for layer_i=1:length(neu_per_layer)-1
            if layer_i==1 && (length(neu_per_layer)-1)>1
                if strcmpi(in_polarity,'unipolar')
                    results_mat=logsig(G_real{layer_i,1}*images_test);
                else
                    results_mat=tansig(G_real{layer_i,1}*images_test);                
                end
            elseif layer_i<length(neu_per_layer)-1 && (length(neu_per_layer)-1)>1
                if strcmpi(in_polarity,'unipolar')
                    results_mat=logsig(G_real{layer_i,1}*results_mat);
                else
                    results_mat=tansig(G_real{layer_i,1}*results_mat);                
                end
            elseif (length(neu_per_layer)-1)>1
                if all_logsig==1
                    if strcmpi(in_polarity,'unipolar')
                        results_mat=logsig(G_real{layer_i,1}*results_mat);
                    else
                        results_mat=tansig(G_real{layer_i,1}*results_mat);
                    end
                else
                    results_mat=G_real{layer_i,1}*results_mat;
                end
                [~,index]=max(results_mat);
                index=index-1;
            else
                results_mat=G_real{layer_i,1}*images_test;
                [~,index]=max(results_mat);
                index=index-1;
            end
        end
        
        comparison=vertcat(labels_test,index);
        comparison(3,:)=comparison(1,:)-comparison(2,:);
        errors=find(comparison(3,:));
        comparison(4,:)=zeros(1,size(comparison,2));
        comparison(4,errors)=1;
        ERROR=sum(comparison(4,:))/size(comparison,2);
        Accuracy_mat=(1-ERROR)*100;
        
        NN_elements.A=A;
        NN_elements.Accuracy=Accuracy;
        NN_elements.CPUtime=tStop;
        %NN_elements.p=p;
        
        if save_pref == 1
            if ~exist(mat_file_weights,'dir')
                mkdir(mat_file_weights)
            end
            save(fullfile(mat_file_weights, 'G_values.mat'),'G','G_real','G_pos','G_neg','G_bin','G_enf','Accuracy','Accuracy_mat','A','labels_test','images_test','tStop');
        end
    end
end

