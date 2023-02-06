function [W_init,W_matrix] = rsCorrect_fcn_mode2(G,mode,digits_order,I_min, I_max,model_ver,memdiode_model,V_read,Rs,partitions)
%RSCORRECT_FCN_MODE2 Summary of this function goes here
%   Detailed explanation goes here

    changes=[];
    prj_dir=fullfile('/home/users','aguirref','nn_rs_uab');

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
        
    G_aux=G;
    G_aux(G_aux<=0)=0;
    G_dual(:,:,1)=G_aux;
    

    G_aux=G;
    G_aux(G_aux>=0)=0;
    G_dual(:,:,2)=-1*G_aux;  
    

    if strcmpi(mode,'log1') || strcmpi(mode,'log2') || strcmpi(mode,'log3') || strcmpi(mode,'log4')

    elseif strcmpi(mode,'lin1') || strcmpi(mode,'lin2') || strcmpi(mode,'lin3') || strcmpi(mode,'lin4')

    elseif strcmpi(mode,'memdiode1') || strcmpi(mode,'memdiode2') || strcmpi(mode,'memdiode3') || strcmpi(mode,'memdiode4')

        if strcmpi(model_ver,'new')
            models_folder=fullfile(prj_dir,'models','memdiode','new_version');
        else
            models_folder=fullfile(prj_dir,'models','memdiode','old_version');
        end

        fid = fopen(fullfile(models_folder,strcat(memdiode_model,'.sp')));
        tline = fgetl(fid);
        while ischar(tline)
%                disp(tline)
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
            if ~isempty(strfind(tline, '.param rsmin=')) || ~isempty(strfind(tline, '+ rsmin='))
                rsmin=str2double(extractAfter(tline,'='));
            end  
            if ~isempty(strfind(tline, '.param rsmax=')) || ~isempty(strfind(tline, '+ rsmax='))
                rsmax=str2double(extractAfter(tline,'='));
            end        
            tline = fgetl(fid);             
        end
        fclose(fid);

        V=V_read;
        repeat_solve=1;
        syms lambda
        partRows=size(G_dual,1)/partitions;

        for i=1:2

            I_matrix_aux = (I_max-I_min)*G_dual(:,:,i) + I_min;        

            I_matrix_aux(I_matrix_aux>I_max)=I_max;
            I_matrix_aux(I_matrix_aux<I_min)=I_min;

            I_matrix(:,:,i)=I_matrix_aux;

            W_matrix(:,:,i) = I_matrix(:,:,i) ;
        end

        while repeat_solve~=0

            for i=1:2

                for G_idx2=1:size(G_dual,2)
                    partRows_i=0;
                    for G_idx1=1:size(G_dual,1)

                        if partRows_i>partRows
                            partRows_i=1;
                        else
                            partRows_i=partRows_i+1;
                        end
                        Rserie=Rs*(partRows_i+G_idx2);
                        eqn=I_matrix(G_idx1,G_idx2,i)==(Imax*lambda+Imin*(1-lambda))*(exp(beta*(alphamax*lambda+alphamin*(1-lambda))*(V-((rsmax*lambda+rsmin*(1-lambda))+Rserie)*I_matrix(G_idx1,G_idx2,i)))-exp(-(1-beta)*(alphamax*lambda+alphamin*(1-lambda))*(V-((rsmax*lambda+rsmin*(1-lambda))+Rserie)*I_matrix(G_idx1,G_idx2,i))))+1e-9;
                        solution=double(vpasolve(eqn,lambda));
                        if ~isempty(solution)
                            W_init(G_idx1,G_idx2,i)=solution;
                            cut_loop=0;
                        else
                            clear W_init
                            I_matrix=I_matrix*0.9;
                            I_matrix(I_matrix>I_max)=I_max;
                            I_matrix(I_matrix<I_min)=I_min;
                            cut_loop=1;
                            break
                        end

                    end
                    if cut_loop==1
                        repeat_solve=1;
                        break;
                    end
                end
                if cut_loop==1
                    repeat_solve=1;
                    break;
                else
                    repeat_solve=0;
                end
            end
        end
    end
end

