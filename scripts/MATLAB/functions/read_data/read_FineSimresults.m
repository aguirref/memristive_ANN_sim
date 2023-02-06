function RAW_DATA = read_FineSimresults(finesim_results)

    if ispc
        current_dir=split(finesim_results,"\");
        current_dir=current_dir(1:end-1);
        result_folder=cell2mat(join(current_dir,"\"));
        result_folder=strrep(root,'\','\\');
    elseif isunix
        current_dir=split(finesim_results,"/");
        current_dir=current_dir(1:end-1);
        result_folder=cell2mat(join(current_dir,"/"));
        %root_str=strrep(root,'\','\\');       
    end

%    parsed_results=strrep(finesim_results,'.tr1','_parsed.tr1');
    
%    fid=fopen(finesim_results);
%     fid_parsed=fopen(parsed_results,'w');
% 
%     in_names=0;
%     in_data=0;
%     line_counter=1;
%     names_line=1e10;
%     data_line=1e10;

%     while ~feof(fid)
%         tline = fgetl(fid);
%         if strfind(tline,'NODES=')
%             nodes_line=strrep(tline,'NODES=','');
%             nodes_line=strrep(nodes_line,'''','');
%             nodes_line=strrep(nodes_line,' ','');
%             nodes_num=str2num(nodes_line);
%             clear nodes_line
%         end
%         if strfind(tline,'#N')
%             in_names=1;
%             in_data=0;
%             if line_counter<names_line
%                 names_line=line_counter;
%             end
%             names=strrep(tline,'#N ','');
% 
%         elseif strfind(tline,'#C')
%             if in_names==1
%                 names=strrep(names,'     ',';');
%                 names=strrep(names,'    ',';');
%                 names=strrep(names,'   ',';'); 
%                 names=strrep(names,'  ',';');        
%                 names=strrep(names,' ',';');  
%                 fprintf(fid_parsed,'%s;\n',names);
%                 line_counter=line_counter+1; 
%                 names_array=split(strrep(names,'''',''),';');
%                 %clear names
%                 in_names=0;
%             end
%             if in_data==1
%                 data=strrep(data,'     ',';');
%                 data=strrep(data,'    ',';');
%                 data=strrep(data,'   ',';'); 
%                 data=strrep(data,'  ',';');        
%                 data=strrep(data,' ',';');              
%                 fprintf(fid_parsed,'%s;\n',data);
%                 line_counter=line_counter+1;            
%                 clear data
%                 in_data=0;
%             end        
%             in_data=1;
%             in_names=0;
%             if line_counter<data_line
%                 data_line=line_counter;
%             end        
%             data=strrep(tline,'#C  ','');
% 
%         else
%             if in_names==1 && in_data==0
%                 names=strcat(names,';',tline);
%             elseif in_names==0 && in_data==1
%                 %fprintf('%s\n',data)                
%                 data=data + "  " + tline;
%             else
%                 fprintf(fid_parsed,'%s\n',tline);
%                 line_counter=line_counter+1;            
%                 %fprintf('%s\n',tline);
%             end
%         end
%     end
%     fclose(fid);
%     %fclose(fid_parsed);

    sim_data=importdata(finesim_results);

    RAW_DATA.variable_mat=sim_data.data(:,2:end)';
    RAW_DATA.time_vect=sim_data.data(:,1)';
    RAW_DATA.variable_name_list=sim_data.colheaders(2:end);
end