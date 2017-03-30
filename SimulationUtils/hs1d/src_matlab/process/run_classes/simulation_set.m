classdef simulation_set
    properties(Access=private)
        mother_folder_directory
        geologic_inputs_directory
        morphologic_inputs_directory
        hydrologic_inputs_directory
        combination_inputs
    end
    
    methods(Access=public)
        function obj=simulation_set(mother_folder_directory)
            obj.mother_folder_directory=mother_folder_directory;
        end
        
        function obj=instantiate_all_inputs_directory(obj)
            obj.geologic_inputs_directory=obj.get_inputs(obj.mother_folder_directory,'GeologicInputs');
            obj.morphologic_inputs_directory=obj.get_inputs(obj.mother_folder_directory,'MorphologicInputs');
            obj.hydrologic_inputs_directory=obj.get_inputs(obj.mother_folder_directory,'HydrologicInputs');
            obj=obj.compute_all_possible_combination_between_inputs;
        end
        
        % find directories in simulations results that encountered errors
        function output_list_directory=get_error_simulation_directory(obj)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            output_list_directory=[];
            compt=1;
            for i=1:length(directory_list)
                A=dir(fullfile(directory,directory_list{i},'*.err'));
                if(~isempty(A))
                    output_list_directory{compt}=fullfile(directory,directory_list{i},A.name);
                    compt=compt+1;
                    if(length(A)>1)
                        fprintf(['WARNING: Potential bad use of .err file chosen: more than one error file found in',fullfile(directory_list{i})]);
                    end
                end
            end
        end
        
        % find directories in simulations results that are incomplete
        function output_list_directory=detect_incomplete_simulations(obj)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            compt=1;
            for i=1:length(directory_list)
                A=dir(fullfile(directory,directory_list{i},'mass_balance.png'));
                if(isempty(A))
                    output_list_directory{compt}=fullfile(directory,directory_list{i});
                    compt=compt+1;
                    if(length(A)>1)
                        fprintf(['WARNING: Potential bad use of file chosen: more than one massbalance file found in',fullfile(directory_list{i})]);
                    end
                end
            end
        end
        
        % find directories in simulation_set that share the same infiltration chronicle
        function output_list_directory=find_simulation_infiltration_chronicle(obj,infiltration_identifier_string)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            length_string=length(infiltration_identifier_string);
            compt=1;
            output_list_directory=[];
            for i=1:length(directory_list)
                Bool=strcmp(directory_list{i}(end-length_string+1:end),infiltration_identifier_string);
                if(Bool==1)
                    output_list_directory{compt}=fullfile(directory,directory_list{i});
                    compt=compt+1;
                end
            end
        end
        
        % find directories in simulation_set that share the same hillslope shape
        function output_list_directory=find_simulation_hillslope_coordinates(obj,Xcoordinates,Ycoordinates)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            compt=1;
            output_list_directory=[];
            for i=1:length(directory_list)
                [X_coordinates_file,Y_coordinates_file]=obj.extract_coordinates(directory_list{i});
                
                Bool1=strcmp(X_coordinates_file,num2str(Xcoordinates));
                Bool2=strcmp(Y_coordinates_file,num2str(Ycoordinates));
                Bool3=Bool1*Bool2;
                if(Bool3==1)
                    output_list_directory{compt}=fullfile(directory,directory_list{i});
                    compt=compt+1;
                end
            end
        end
        
        function output_list_directory=find_simulation_hillslope_sloptype(obj,sloptype_str)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            compt=1;
            output_list_directory=[];
            for i=1:length(directory_list)
                [sloptype]=obj.extract_slope_type(directory_list{i});
                
                Bool=strcmp(sloptype,sloptype_str);
                if(Bool==1)
                    output_list_directory{compt}=fullfile(directory,directory_list{i});
                    compt=compt+1;
                end
            end
        end
        
        function [X_coordinates_file,Y_coordinates_file]=extract_coordinates(obj,directory_path_name)
            Id1=strfind(directory_path_name,'X');
            Id2=strfind(directory_path_name,'Y');
            Id3=strfind(directory_path_name,'slop');
            X_coordinates_file=directory_path_name(Id1+2:Id2-2);
            Y_coordinates_file=directory_path_name(Id2+2:Id3-2);
        end
        
        function [sloptype]=extract_slope_type(obj,directory_path_name)
            Id3=strfind(directory_path_name,'slop');
            sloptype=directory_path_name(Id3+4:Id3+6);
        end
                
        % find directories in simulation_set that share the same hillslope shape
        function output_list_directory=find_simulation_geometric_properties(obj,f,k,d)
            directory_list=obj.get_output_simulation;
            directory=[obj.mother_folder_directory,'\Simulations'];
            compt=1;
            output_list_directory=[];
            for i=1:length(directory_list)
                directory_total=fullfile(directory,directory_list{i},'geologic.input');
                fid=fopen(directory_total,'r');
                if(fid>0)
                    C = textscan(fid, '%s','delimiter', '\t');
                    f_file=str2num(C{1}{end-2});
                    k_file=str2num(C{1}{end-1});
                    d_file=str2num(C{1}{end});
                    Bool=(f_file-f)^2+(k_file-k)^2+(d_file-d)^2;
                    
                    if(Bool==0)
                        output_list_directory{compt}=fullfile(directory,directory_list{i});
                        compt=compt+1;
                    end
                    fclose(fid);
                end
            end
        end
    end
    
    methods(Access=private)    
        function output_list_directory=get_inputs(obj,directory,type_of_inputs)
            if(nargin>=3)
                directory=[directory,'\',type_of_inputs];
            end
            list=list_folder_of(directory);
            output_list_directory=[];
            compt=1;
            for i=1:length(list)
                A=dir(fullfile(directory,list{i},'*.input'));
                if(~isempty(A))
                    output_list_directory{compt}=fullfile(directory,list{i},A.name);
                    compt=compt+1;
                    if(length(A)>1)
                        fprintf(['WARNING: Potential bad use of .input file chosen: more than one file in',fullfile(directory,list{i})]);
                    end
                end
            end
        end
        
        function simulation_list_directory=get_output_simulation(obj)
            directory=obj.mother_folder_directory;
            directory=[directory,'\Simulations'];
            simulation_list_directory=list_folder_of(directory);
        end
        
        function obj=compute_all_possible_combination_between_inputs(obj)
            obj.combination_inputs=allcomb(obj.geologic_inputs_directory,obj.hydrologic_inputs_directory,obj.morphologic_inputs_directory);
        end
        
        function obj=run_simulation_set(obj)
            % create a Simulation folder in the mother_folder_directory
            simulation_folder_root=fullfile(obj.mother_folder_directory,'Simulations');
            folder_create(simulation_folder_root);
%             obj.initialize_summary_file(simulation_folder_root);
            % launch in a for loop all the simulations
            max_iter=size(obj.combination_inputs);
            max_iter=max_iter(1);
            
%             perm_pos=randperm(max_iter);
            
            % generate some pdfs with wich the runs will update f k d
            % ##JM maybe to change - to indroduce in new versions
%             pd_f=makedist('Normal','mu',0.3,'sigma',0.1);
%             t_f=truncate(pd_f,0.05,0.5);
%             pd_k=makedist('Lognormal','mu',1,'sigma',1);
%             t_k=truncate(pd_k,0.05,15);
%             pd_d=makedist('LogNormal','mu',1,'sigma',1);
%             t_d=truncate(pd_d,0.2,11);
%             r_f=random(t_f,max_iter,1);
%             r_k=random(t_k,max_iter,1);
%             r_d=random(t_d,max_iter,1);
            
            for i=1:max_iter               
                % locations of different inputs file
                geo_loc=obj.combination_inputs(i,1);
                hydro_loc=obj.combination_inputs(i,2);
                morpho_loc=obj.combination_inputs(i,3);
                 % ##JM maybe to change - to indroduce in new versions
%                 geo_loc=obj.combination_inputs(perm_pos(i),1);
%                 hydro_loc=obj.combination_inputs(perm_pos(i),2);
%                 morpho_loc=obj.combination_inputs(perm_pos(i),3);
%                 AA=obj.combination_inputs(i,3);
%                 M2=obj.read_input_file(AA{1});
%                 w_test=M2(:,2);
%                 if(w_test(1)<25)
                % 1/ create a specific folder to store parameters and results files of the run
                %    copy paste the inputs file and some figures in the simulation results folder
%                 for j=1:61
                [destination_geol_file,destination_hydro_file,destination_morpho_file,folder_output]=obj.create_specific_results_folder(simulation_folder_root,morpho_loc,hydro_loc,geo_loc);

                
                % 2/ read input files and "recreate" objects (hs1D and source objects)
                % ##JM maybe to change - to indroduce in new versions
                M=obj.read_input_file(destination_geol_file);
                %#JM change after
                f=M(1,1); k=M(1,2); d=M(1,3);
                % ##JM maybe to change - to indroduce in new versions
%                 f=r_f(i); k=r_k(i); d=r_d(i);
%                 dlmwrite(destination_geol_file,[f,k,d],'-append','delimiter','\t','precision','%E');
                M=obj.read_input_file(destination_morpho_file);
                %#JM change after test
                x=M(:,1); w=M(:,2); slope_angle=M(:,3);
%                 if(j==31)
%                     slope_angle=0.15*ones(size(x));
%                 else
%                     slope_angle=((((j-1)/2):(62-2*j)/(2*(length(x)-1)):((61-j)/2))/100)';
%                 end
%                 slope_angle(x>=100)=0.05;
%                 slope_angle(x>150)=0.1;
%                 slope_angle=0.01*j*ones(size(slope_angle));
%                     w(w<25)=30;
                    % ##JM maybe to change - to indroduce in new versions
%                 hs1D=hillslope1D(perm_pos(i),f,k);
                hs1D=hillslope1D(i,f,k);
%                 hs1D=hs1D.set_spatial_parameters(x,w,slope_angle,d*ones(size(x)));
%                 d=(0.5:1/(length(x)-1):1.5)';
%                 d=ones(size(x));
%                 d=0.5*j*ones(size(x));
%                 d(x<50)=0.5;
%                 d=0.5*ones(size(x));
                hs1D=hs1D.set_spatial_parameters(x,w,slope_angle,d);
                hs1D=hs1D.compute_elevation;
                % #JM maybe to change after
%                 hs1D.plot_save_width_function(folder_output);
                hs1D.plot_save_elevation_function(folder_output);
                hs1D.plot_save_slope_angle_function(folder_output);
                [M,input_type]=obj.read_input_file(destination_hydro_file);
                t=M(:,1);
                recharge_chronicle=(M(:,2))';
                ratio_P_R=1;
                source_terms=source('data_based');
                [~,source_terms]=source_terms.set_recharge_chronicle_data_based(t/(3600*24),ratio_P_R,recharge_chronicle,'m/s');
                
                % 3/ create a runs object and run the simulation
                run_obj=runs;
                % set the solver options default or assigned in parameters via an odeset structure
                % specify Refine options for real infiltrations chronicle because for accuracy you need 
                % to force matlab ode15s to compute where you know sthg is happening
                if(strcmp(input_type,'Real1') || strcmp(input_type,'Real2') || strcmp(input_type,'Real3'))
%                     odeset_struct=odeset('RelTol',2.5e-14,'AbsTol',1e-14,'Refine',-1);
                    odeset_struct=odeset('RelTol',2.5e-14,'Refine',-1);
                % otherwise let matlab define its own time range for efficiency reasons (speed)
                else
%                     odeset_struct=odeset('RelTol',2.5e-14,'AbsTol',1e-14);
                    odeset_struct=odeset('RelTol',2.5e-14);
                end
                solver_options=run_obj.set_solver_options(odeset_struct);

                % if the simulation is not already a steady state in itself...
                if(~strcmp(input_type,'Synthetic1'))
                    % ... preprocess simulations to reach a steady state with averaged forcing conditions
                    % ##JM better change after
%                     recharge_averaged=0;
                    recharge_averaged=1e3*24*3600*source_terms.recharge_mean; % recharge averaged in mm/d
                    state_values_initial=obj.prerun_steady_state(hs1D,recharge_averaged);
                    presteadystate_percentage_loaded=-2; % -2 is the key to start a simulation with a customed initial condition for storage prescribed in Sinitial
                    % run the simulation starting from the steady state condition 
                    run_obj=run_obj.run_simulation(hs1D,source_terms,presteadystate_percentage_loaded,solver_options,state_values_initial);
                else
                    % run the simulation starting from empty hillslope
                    percentage_loaded=0;
                    run_obj=run_obj.run_simulation(hs1D,source_terms,percentage_loaded,solver_options);
                end
                tmax=t(end);
                error=run_obj.boussinesq_simulation.save_error_file(tmax,folder_output);
                run_obj.boussinesq_simulation.save_info_file(folder_output);
                
                % save boussinesq_simulation & hs1D objects
                space_limited_code=2; % if no code, assumed that enough space to store directly .mat objects, if code=1, enough space to convert it in .txt, if code=2, ode15s details are not saved
                run_obj.save_sim_run(folder_output,space_limited_code);
                
                % save simulation_results objects and simulation_results output files
                t_sim_results=run_obj.simulation_results.t;
                
                [S_max,Watershed_area,w,Flux_in_total,Watershed_area_spatialized]=run_obj.save_key_parameters(t_sim_results,folder_output);
                run_obj.simulation_results.save_results(S_max,Watershed_area,folder_output);
                run_obj.simulation_results.plot_results(S_max,Watershed_area,w,Flux_in_total,folder_output);
%                 run_obj.simulation_results.plot_mass_changes(Watershed_area,Flux_in_total',Watershed_area_spatialized,folder_output);
                dS=run_obj.boussinesq_simulation.compute_raw_storage_variation;
                Flux_in=run_obj.compute_flux_in(t_sim_results);
                size_dS=size(dS);
                dS_rough=nan(size_dS(1),length(t_sim_results));
                for k=1:size_dS(1)
                    dS_rough(k,:)=interpn(run_obj.boussinesq_simulation.sol_simulated.x,dS(k,:),t_sim_results);
                end
                [Mass_balance_tot_mmd,Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,Storage_Var_mmd]= ...
                    run_obj.simulation_results.mass_changes_study(Watershed_area_spatialized,Flux_in,dS_rough);
                run_obj.simulation_results.plot_mass_balance(Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,Storage_Var_mmd,Mass_balance_tot_mmd,folder_output);
                 % ##JM maybe to change - to indroduce in new versions
%                 obj.append_summary_file(simulation_folder_root,perm_pos(i),folder_output,error);
                obj.append_summary_file(simulation_folder_root,i,folder_output,error);
                clear run_obj;
%                 end
            end
        end
        
        function state_values_initial=prerun_steady_state(obj,hs1D,recharge_averaged)
            percentage_loaded=0;
            source_steady=source('steady');
            tsteady=0:1:3650;
            tmin=tsteady(1); tmax=tsteady(end); Nt=length(tsteady); time_unity_type='day';
            t=time_properties(tmin,tmax,Nt,time_unity_type);
            source_steady=source_steady.set_recharge_chronicle(t,recharge_averaged);
            prerun_steady=runs;
            % add an event equal to 1 to detect steady_state
            detect_steady_state=1;
%             odeset_struct=odeset('RelTol',2.5e-14,'AbsTol',1e-17,'Events',detect_steady_state);
            odeset_struct=odeset('RelTol',2.5e-14,'Events',detect_steady_state);
            solver_options=prerun_steady.set_solver_options(odeset_struct);
            prerun_steady=prerun_steady.run_simulation(hs1D,source_steady,percentage_loaded,solver_options);
%             Sinitial=prerun_steady.get_final_storage;
            state_values_initial=prerun_steady.get_final_state_values;
        end
        
        function [M,input_type]=read_input_file(obj,filename)
            fid = fopen(filename, 'r');
            if(fid>0)
                input_type= textscan(fid,'%s', 1,'headerlines',1);
                input_type=input_type{1};
                M=dlmread(filename,'\t',3,0);
                fclose(fid);
            else
                fprintf('');
            end
        end
        
        function copypaste_png_figures(obj,file_location,folder_sep_pos,folder_output)
            png_directories=dir(fullfile(file_location{1}(1:folder_sep_pos{1}(end)-1),'*.png'));
            for k=1:length(png_directories)
                filename_png=fullfile(file_location{1}(1:folder_sep_pos{1}(end)-1),png_directories(k).name);
                copypaste_file_customed(filename_png,folder_output);
            end
        end
        
        function [destination_geol_file,destination_hydro_file,destination_morpho_file,folder_output]=create_specific_results_folder(obj,simulation_folder_root,morpho_loc,hydro_loc,geo_loc)
            % 1/ create a specific folder to store parameters and results files of the run
            c=clock; time_string_folder=strcat(num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)),'_',num2str(round(c(6))));
            folder_sep_pos=strfind(morpho_loc,'\');
            coordinate_string=morpho_loc{1}(folder_sep_pos{1}(end-1)+1:folder_sep_pos{1}(end)-1);
            folder_sep_pos_hydro=strfind(hydro_loc,'\');
            hydro_type_string=hydro_loc{1}(folder_sep_pos_hydro{1}(end-1)+1:folder_sep_pos_hydro{1}(end)-1);
            folder_sep_pos_geo=strfind(geo_loc,'\');
            geo_type_string=geo_loc{1}(folder_sep_pos_geo{1}(end-1)+1:folder_sep_pos_geo{1}(end)-1);
            % only the hydro and the morpho identifiers in folders name
            name_folder=[time_string_folder,'_',coordinate_string,'_',hydro_type_string];
            folder_output=fullfile(simulation_folder_root,name_folder);
            folder_create(folder_output);
            
            % 2/ write an parameters file with the different input properties
            obj.write_input_file_parameters(folder_output,coordinate_string,hydro_type_string,geo_type_string);
            
            % 3/ copy paste the inputs file and some figures and the simulation results folder
            % get png images in simulations folder (previously in morphologic & hydrologic inputs directory)
            obj.copypaste_png_figures(morpho_loc,folder_sep_pos,folder_output);
            obj.copypaste_png_figures(hydro_loc,folder_sep_pos_hydro,folder_output);
            % copy paste the parameter folders in the simulation file folder file
            destination_geol_file=copypaste_file_customed(geo_loc{1},folder_output);
            destination_hydro_file=copypaste_file_customed(hydro_loc{1},folder_output);
            destination_morpho_file=copypaste_file_customed(morpho_loc{1},folder_output);
        end
        
        function write_input_file_parameters(obj,foldername,coordinate_string,hydro_type_string,geo_type_string)
            folder_sep_pos_coord=strfind(coordinate_string,'_');
            X_coordinates=coordinate_string(folder_sep_pos_coord(1)+1:folder_sep_pos_coord(2)-1);
            Y_coordinates=coordinate_string(folder_sep_pos_coord(3)+1:folder_sep_pos_coord(4)-1);
            folder_sep_pos=strfind(geo_type_string,'_');
            f_values=geo_type_string(folder_sep_pos(1)+1:folder_sep_pos(2)-1);
            k_values=geo_type_string(folder_sep_pos(3)+1:folder_sep_pos(4)-1);
            d_values=geo_type_string(folder_sep_pos(5)+1:end);

            filename_input=strcat(foldername,'\input.param');
            fid = fopen(filename_input, 'w');
            
            fprintf(fid,'Hydrologic Forcing Type \n');
            fprintf(fid,[hydro_type_string,'\n']);
            fprintf(fid,'X\tY\n');
            fprintf(fid,[X_coordinates,'\t',Y_coordinates,'\n']);
            fprintf(fid,'f\tk\td\n');
            fprintf(fid,[f_values,'\t',k_values,'\t',d_values,'\n']);
            fclose(fid);
        end
        
        function initialize_summary_file(obj,simulation_folder_root)
            filename=strcat(simulation_folder_root,'\summary_file.txt');
            fid = fopen(filename, 'w');
            fprintf(fid, 'SimId\tFolderPath\tErr\n');
            fprintf(fid, '\n');
            fclose(fid);
        end
        
        function append_summary_file(obj,simulation_folder_root,i,simulation_folder_path,error)
            filename=strcat(simulation_folder_root,'\summary_file.txt');
            fid = fopen(filename, 'a');
            simulation_folder_path=strrep(simulation_folder_path, '\', '/');
            fprintf(fid, [num2str(i),'\t',simulation_folder_path,'\t',num2str(error),'\n']);
            fclose(fid);
        end
    end
    
    methods(Static)
        function run_simulations(folder_root)
            sim_set=simulation_set(folder_root);
            sim_set=sim_set.instantiate_all_inputs_directory;
            sim_set=sim_set.run_simulation_set;
        end
        
        function rerun_simulation(SimulationPath)
            sim_set=simulation_set(SimulationPath);
            output_list_directory=dir(fullfile(SimulationPath,'*.input'));
            sim_set.geologic_inputs_directory{1}=fullfile(SimulationPath,output_list_directory(1).name);
            sim_set.morphologic_inputs_directory{1}=fullfile(SimulationPath,output_list_directory(3).name);
            sim_set.hydrologic_inputs_directory{1}=fullfile(SimulationPath,output_list_directory(2).name);
            sim_set=sim_set.compute_all_possible_combination_between_inputs;
            sim_set.run_simulation_set;
        end
        
        function run_simulations2(hs1D,folder_root,n)
            % generate ad hoc parametrization defined in geologic_input_set static methods
            [geol_input_set,val]=geologic_input_set.generate_customed_hillslope_parametrization;
            % generate sistematic parametrization defined in generate_systematic_hillslope_parametrization methods
            %             [geol_input_set,val]=geologic_input_set.generate_systematic_hillslope_parametrization(n);
            folder_root=strcat(folder_root,'\Simulations\');
            for i=1:length(val)
                % create folder output
                c=clock; time_string_folder=strcat(num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)),'_',num2str(c(6)));
                folder_output=strcat(folder_root,time_string_folder);
                folder_create(folder_output);
                simulation_run=run_set;
                [simulation_run,discretization,percentage_loaded]=simulation_run.choose_structure_parameters(hs1D);
                if(size(val,2)==3)
                    simulation_run=simulation_run.choose_hydraulic_parameters(val(i,1),val(i,2),val(i,3));
                elseif(size(val,2)==4)
                    simulation_run=simulation_run.choose_hydraulic_parameters(val(i,1),val(i,2),val(i,3),val(i,4));
                end
                [t,source_term]=simulation_run.set_source_options;
                boundary_cond=simulation_run.set_boundary_conditions(discretization);
                solver_options=simulation_run.set_solver_options;
                simulation_run=simulation_run.initialize_boussinesq_simulation(discretization,source_term,boundary_cond,percentage_loaded,t.get_tmin);
                simulation_run.save_parametrization(folder_output);
                % Run the DAE
                [t,x_center,x_edge,S,Q,QS,obj.boussinesq_simulation]=simulation_run.boussinesq_simulation.implicit_scheme_solver(t,solver_options);
                simulation_run.simulation_results=simulation_results(t,x_center,x_edge,S,Q,QS);
                simulation_run.save_results(folder_output);
            end
        end
    end
end