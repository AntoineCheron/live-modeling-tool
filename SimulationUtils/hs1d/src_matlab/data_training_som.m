classdef data_training_som
    properties(Access=public)
        max_length
        width_disc
        elev_disc
        hydraulic_param
        matrix_train
        
        Area_param
        angle_mean_param
        width_origin_param
        
        Length_sat_mean
        Length_sat_std
        Area_sat_mean
        Area_sat_std
    end
    
    methods(Access=public)
        function obj=data_training_som(max_length,width_disc,elev_disc,hydraulic_param,Area_param,angle_mean_param,width_origin_param,Length_sat_mean,Length_sat_std,Area_sat_mean,Area_sat_std)
            obj.max_length=max_length;
            obj.width_disc=width_disc;
            obj.elev_disc=elev_disc;
            obj.hydraulic_param=hydraulic_param;
%             size_max_length=size(max_length);
%             size_width_disc=size(width_disc);
%             size_elev_disc=size(elev_disc);
            size_hydraulic_param=size(hydraulic_param);
%             test_null=(size_max_length(1)-size_width_disc(1))^2+(size_max_length(1)-size_elev_disc(1))^2+(size_max_length(1)-size_hydraulic_param(1))^2;
%             max_length_normalized=max_length/max(max_length);
%             width_disc_normalized=bsxfun(@rdivide,width_disc,max(width_disc));
%             elev_disc_normalized=bsxfun(@rdivide,elev_disc,max(elev_disc));
%             hydraulic_param_normalized=bsxfun(@rdivide,hydraulic_param,max(hydraulic_param));
%             if(test_null==0)
%                 obj.matrix_train=nan(size_max_length(1),size_max_length(2)+size_width_disc(2)+size_elev_disc(2)+size_hydraulic_param(2));
%                 obj.matrix_train(1:size_max_length(1),1:size_max_length(2))=max_length_normalized;
%                 obj.matrix_train(1:size_max_length(1),size_max_length(2)+1:size_max_length(2)+size_width_disc(2))=width_disc_normalized;
%                 obj.matrix_train(1:size_max_length(1),size_max_length(2)+size_width_disc(2)+1:size_max_length(2)+size_width_disc(2)+size_elev_disc(2))=elev_disc_normalized;
%                 obj.matrix_train(1:size_max_length(1),size_max_length(2)+size_width_disc(2)+size_elev_disc(2)+1:size_max_length(2)+size_width_disc(2)+size_elev_disc(2)+size_hydraulic_param(2))=hydraulic_param_normalized;
%             end
%             
                        
            obj.Area_param=Area_param;
            obj.angle_mean_param=angle_mean_param;
            obj.width_origin_param=width_origin_param;
            
            obj.Length_sat_mean=Length_sat_mean;
            obj.Length_sat_std=Length_sat_std;
            obj.Area_sat_mean=Area_sat_mean;
            obj.Area_sat_std=Area_sat_std;
            
            saturated_param=[Area_sat_mean,Area_sat_std];
            size_saturated_param=size(saturated_param);
            size_Area=size(Area_param);
            size_angle_mean=size(angle_mean_param);
            size_width_origin=size(width_origin_param);
            test_null=(size_Area(1)-size_width_origin(1))^2+(size_Area(1)-size_angle_mean(1))^2+(size_Area(1)-size_hydraulic_param(1))^2;
            if(test_null==0)
                obj.matrix_train=nan(size_width_origin(1),size_width_origin(2)+size_Area(2)+size_angle_mean(2)+size_hydraulic_param(2));
                obj.matrix_train(1:size_width_origin(1),1:size_width_origin(2))=width_origin_param;
                obj.matrix_train(1:size_width_origin(1),size_width_origin(2)+1:size_width_origin(2)+size_Area(2))=Area_param;
                obj.matrix_train(1:size_width_origin(1),size_width_origin(2)+size_Area(2)+1:size_width_origin(2)+size_Area(2)+size_angle_mean(2))=angle_mean_param;
                obj.matrix_train(1:size_width_origin(1),size_width_origin(2)+size_Area(2)+size_angle_mean(2)+1:size_width_origin(2)+size_Area(2)+size_angle_mean(2)+size_hydraulic_param(2))=hydraulic_param;
                obj.matrix_train(1:size_width_origin(1),size_width_origin(2)+size_Area(2)+size_angle_mean(2)+size_hydraulic_param(2)+1:size_width_origin(2)+size_Area(2)+size_angle_mean(2)+size_hydraulic_param(2)+size_saturated_param(2))=saturated_param;
                
            end
        end
        
        function obj=normalize_matrix(obj)
            mean_matrix=mean(obj.matrix_train);
            std_matrix=std(obj.matrix_train);
            obj.matrix_train=bsxfun(@minus,obj.matrix_train,mean_matrix);
            obj.matrix_train=bsxfun(@rdivide,obj.matrix_train,std_matrix);
            
        end
    end
end