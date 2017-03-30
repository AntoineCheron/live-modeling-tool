classdef geologic_input_set
    properties(Access=public)
        n
        f % {[fmin, fmax], 'type'} drainable porosity
        k % {[kmin, kmax], 'type'} hydraulic 
        d % {[dmin, dmax], 'type'} soil_depth
        i % optional {[imin, imax], 'type'} slope angle
        parameters_grid_list
    end
    
    methods(Static)
        
        function [geol_input_set,val]=generate_systematic_hillslope_parametrization(n)
            f={[0.1,0.5],'lin'};
            k={[1e-5,1],'log'};
            d={[0.5,2],'lin'};
            
            geol_input_set=geologic_input_set(f,k,d,n);
%             i={[0.01,0.3],'lin'};
%             geol_input_set=geologic_input_set(f,k,d,n,i);
            [geol_input_set,val]=geol_input_set.set_different_parametrization;
        end
        
        function [geol_input_set,val]=generate_customed_hillslope_parametrization(file_directory)
            f={[0.05,0.1,0.3,0.4],'custom'};
            k={[0.001,0.01,0.1,1,10],'custom'};
            d={[0.5,1,2,10],'custom'};
            nparam=nan;
            geol_input_set=geologic_input_set(f,k,d,nparam);
            [geol_input_set,val]=geol_input_set.set_different_parametrization;
            if(nargin>0)
                geol_input_set.save_parametrization(val,file_directory);
                geol_input_set.save_mat_file(file_directory);
            end
        end
    end
    
    methods(Access=private)
        function obj=geologic_input_set(f,k,d,nparam,i)
            if(nargin<5) 
                obj.i=[];
            else obj.i=i;
            end
            obj.n=nparam;
            obj.f=f;
            obj.k=k;
            obj.d=d;
        end
        
        function [obj,val]=set_different_parametrization(obj)
            if(~isnan(obj.n))
                n_param=round((obj.n)^(1/4));
            else
                n_param=-1;
            end
            p{1}=parameter_grid(obj.f{1}(1),obj.f{1}(2:end),n_param,strcmp(obj.f{2},'log')-strcmp(obj.f{2},'custom'));
            p{2}=parameter_grid(obj.k{1}(1),obj.k{1}(2:end),n_param,strcmp(obj.k{2},'log')-strcmp(obj.k{2},'custom'));
            p{3}=parameter_grid(obj.d{1}(1),obj.d{1}(2:end),n_param,strcmp(obj.d{2},'log')-strcmp(obj.d{2},'custom'));
            if(~isempty(obj.i))
                p{4}=parameter_grid(obj.i{1}(1),obj.i{1}(2:end),n_param,strcmp(obj.i{2},'log')-strcmp(obj.i{2},'custom'));
            end
            obj.parameters_grid_list=parameter_grid_list(p);
            val=obj.parameters_grid_list.cross_product;
        end
        
        function save_parametrization(obj,val,file_directory)
            for i=1:length(val)
                if(size(val,2)==4)
                    folder_name_string=['f_',num2str(val(i,1)),'_k_',num2str(val(i,2)),'_d_',num2str(val(i,3)),'_i_',num2str(val(i,4))];
                    string_char=sprintf('f[percent]\tk[m/h]\td[m]\ti[percent]\n');
                    string_char2=sprintf('Geologic data include slope angle i. f is the drainable porosity. k is the hydraulic conductivity. d is soil depth. i is slope angle. \n');
                elseif(size(val,2)==3)
                    folder_name_string=['f_',num2str(val(i,1)),'_k_',num2str(val(i,2)),'_d_',num2str(val(i,3))];
                    string_char=sprintf('f[percent]\tk[m/h]\td[m]\n');
                    string_char2=sprintf('Geologic data does not include slope angle i. f is the drainable porosity. k is the hydraulic conductivity. d is soil depth. \n');
                else
                    fprintf('WARNING: Potential problem with the size of val');
                end
                file_directory_final=fullfile(file_directory,'GeologicInputs',folder_name_string);
                folder_create(file_directory_final);
                filename=[file_directory_final,'\geologic.input'];
                M=val(i,:);
                fid = fopen(filename, 'w');
                string_char3=sprintf('Customed geologic data taken \n');
                fprintf(fid, string_char3);
                fprintf(fid, string_char2);
                fprintf(fid, string_char);
                fclose(fid);
                dlmwrite(filename,M, '-append', 'precision', '%E','delimiter','\t');
                save(fullfile(file_directory_final,'geologic_parametrization_matrix.mat'),'val');
            end
        end
        
        function save_mat_file(obj,file_directory)
            save(fullfile(file_directory,'GeologicInputs','geologic_input_set.mat'),'obj');
        end
    end
end