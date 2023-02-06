function [W_init,W_matrix] = G2W_map(G,H_max,H_min,mode,digits_order,I_min, I_max,model_ver,memdiode_model,V_read,variability, write, varargin)
    changes=[];
    prj_dir=fullfile('/home/users','aguirref','nn_rs_uab');
    layer_i=0;
    total_layers=0;
    
    %% Optional values assignment
    if mod(nargin-12,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-12
            if strcmp(varargin{i},'layer_i')
                layer_i=varargin{i+1};
            end
            if strcmp(varargin{i},'total_layers')
                total_layers=varargin{i+1};
            end
            if strcmp(varargin{i},'prj_dir')
                prj_dir=varargin{i+1};
            end  
            if strcmp(varargin{i},'ANN_type')
                ANN_type=varargin{i+1};
            end  
            if strcmp(varargin{i},'simEngine')
                simEngine=varargin{i+1};
            end
        end
    end

    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
        model_syntax_ver='HSPICE';
    elseif strcmp(simEngine,'LTspice')
        model_syntax_ver='LTspice';
    else
        model_syntax_ver='HSPICE';
    end
    
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
    
    if size(G,3)==1
        G_aux=G';
        G_aux(G_aux<=0)=0;
        G_dual(:,:,1)=G_aux;


        G_aux=G';
        G_aux(G_aux>=0)=0;
        G_dual(:,:,2)=-1*G_aux;  
    else
        G_dual(:,:,1)=squeeze(G(1,:,:))';
        G_dual(:,:,2)=squeeze(G(2,:,:))';
    end
    
    if H_max<=H_min
        aux=H_max;
        H_max=H_min;
        H_min=aux;
        fprintf('the input values for H_max and H_min are not correct: H_min should be smaller than H_max. By default they have been switched over');
    end
    
    solved_eq=0;    
    if total_layers>0
        label_str=sprintf('Calculating CPA weights for layer %d', layer_i);
    end
    
    for i=1:2
        if strcmpi(ANN_type,'BNN')
            W_init_aux(:,:) = G_dual(:,:,i)+normrnd(0,variability,[size(G_dual,1) size(G_dual,2)]);
            W_init_aux(W_init_aux>1)=1;
            W_init_aux(W_init_aux<0)=0;
            W_init(:,:,i)=W_init_aux;
            W_matrix(:,:,i) = (I_max-I_min)*G_dual(:,:,i) + I_min;
        else
            if strcmpi(mode,'log1') || strcmpi(mode,'log2') || strcmpi(mode,'log3') || strcmpi(mode,'log4')
                W_init(:,:,i) = 10.^((log10(H_max)-log10(H_min))*G_dual(:,:,i) + log10(H_min));   
            
            elseif strcmpi(mode,'lin1') || strcmpi(mode,'lin2') || strcmpi(mode,'lin3') || strcmpi(mode,'lin4')
                W_init(:,:,i) = (H_max-H_min)*G_dual(:,:,i) + H_min;   
                W_matrix(:,:,i) = (I_max-I_min)*G_dual(:,:,i) + I_min;
            
            elseif strcmpi(mode,'memdiode1') || strcmpi(mode,'memdiode2') || strcmpi(mode,'memdiode3') || strcmpi(mode,'memdiode4')
                
                if strcmpi(model_ver,'DMM')
                    models_folder=fullfile(prj_dir,'models','memdiode','DMM',model_syntax_ver);
                elseif strcmpi(model_ver,'QMM')
                    models_folder=fullfile(prj_dir,'models','memdiode','QMM',model_syntax_ver);
                elseif strcmpi(model_ver,'QMM_SBSF')
                    models_folder=fullfile(prj_dir,'models','memdiode','QMM_SBSF',model_syntax_ver);
                end
                
                fid = fopen(fullfile(models_folder,strcat(memdiode_model,'.sp')));
                tline = fgetl(fid);
                while ischar(tline)
    %                disp(tline)
                    if ~isempty(strfind(tline, '.param beta=')) || ~isempty(strfind(tline, '+ beta='))
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
                    if ~isempty(strfind(tline, '.param rsmin=')) || ~isempty(strfind(tline, '+ rsmin='))
                        rsmin=str2double(extractAfter(tline,'='));
                    end  
                    if ~isempty(strfind(tline, '.param rsmax=')) || ~isempty(strfind(tline, '+ rsmax='))
                        rsmax=str2double(extractAfter(tline,'='));
                    end        
                    if ~isempty(strfind(tline, '.param imax_intercept=')) || ~isempty(strfind(tline, '+ imax_intercept='))
                        imax_intercept=str2double(extractAfter(tline,'='));
                    end  
                    if ~isempty(strfind(tline, '.param imax_slope=')) || ~isempty(strfind(tline, '+ imax_slope='))
                        imax_slope=str2double(extractAfter(tline,'='));
                    end        
                    if ~isempty(strfind(tline, '.param imin_intercept=')) || ~isempty(strfind(tline, '+ imin_intercept='))
                        imin_intercept=str2double(extractAfter(tline,'='));
                    end  
                    if ~isempty(strfind(tline, '.param imin_slope=')) || ~isempty(strfind(tline, '+ imin_slope='))
                        imin_slope=str2double(extractAfter(tline,'='));
                    end  
                    if ~isempty(strfind(tline, '.param alphamax_intercept=')) || ~isempty(strfind(tline, '+ alphamax_intercept='))
                        alphamax_intercept=str2double(extractAfter(tline,'='));
                    end  
                    if ~isempty(strfind(tline, '.param alphamax_slope=')) || ~isempty(strfind(tline, '+ alphamax_slope='))
                        alphamax_slope=str2double(extractAfter(tline,'='));
                    end        
                    if ~isempty(strfind(tline, '.param alphamin_intercept=')) || ~isempty(strfind(tline, '+ alphamin_intercept='))
                        alphamin_intercept=str2double(extractAfter(tline,'='));
                    end  
                    if ~isempty(strfind(tline, '.param alphamin_slope=')) || ~isempty(strfind(tline, '+ alphamin_slope='))
                        alphamin_slope=str2double(extractAfter(tline,'='));
                    end                      
                    if ~isempty(strfind(tline, 'x0='))
                        oxygen_flow=str2double(extractAfter(tline,'x0='));
                    end  
                    tline = fgetl(fid);             
                end
                fclose(fid);
                
                if isnan(Imax)
                    Imax=imax_slope*oxygen_flow+imax_intercept;
                end
                if isnan(Imin)
                    Imin=imin_slope*oxygen_flow+imin_intercept;
                end
                if isnan(alphamax)
                    alphamax=alphamax_slope*oxygen_flow+alphamax_intercept;
                end
                if isnan(alphamin)
                    alphamin=alphamin_slope*oxygen_flow+alphamin_intercept;
                end

                V=V_read;
                
                I_matrix_aux = (I_max*0.9-I_min)*G_dual(:,:,i) + I_min;        
    
                I_matrix_aux = I_matrix_aux.*normrnd(1,variability,[size(G_dual,1) size(G_dual,2)]);  
                
                I_matrix_aux(I_matrix_aux>I_max*0.9)=I_max*0.9;
                I_matrix_aux(I_matrix_aux<I_min)=I_min;
                
                I_matrix(:,:,i) = I_matrix_aux;
                
                W_matrix(:,:,i) = I_matrix(:,:,i);
                
                syms lambda
                
                if write==0
                    for G_idx2=1:size(G_dual,2)
                        for G_idx1=1:size(G_dual,1)
    
                            eqn=I_matrix(G_idx1,G_idx2,i)==(Imax*lambda+Imin*(1-lambda))*(exp(beta*(alphamax*lambda+alphamin*(1-lambda))*(V-((rsmax*lambda+rsmin*(1-lambda)))*I_matrix(G_idx1,G_idx2,i)))-exp(-(1-beta)*(alphamax*lambda+alphamin*(1-lambda))*(V-((rsmax*lambda+rsmin*(1-lambda)))*I_matrix(G_idx1,G_idx2,i))));
                            W_init(G_idx1,G_idx2,i)=double(vpasolve(eqn,lambda));   
                            if total_layers>0
                                solved_eq=solved_eq+1;
                                multiWaitbar(label_str,solved_eq/(size(G_dual,1)*size(G_dual,2)*2));
                            end
                        end    
                    end
                else
                    W_init(:,:,i) = zeros(size(I_matrix(:,:,i)));
                end
            end
        end
    end
    fh = findobj( 'Type', 'Figure', 'Name', 'memdiode_IV_characteristic' );
    close(fh);
end

