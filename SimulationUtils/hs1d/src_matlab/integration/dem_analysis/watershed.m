classdef watershed
    % Set of hillslope
    properties(Access=public)
        DEM            % GRIDobj topotoolbox
        Outlet         % x y coordinates and ID (from topotoolbox) (1x3 array)
        Confluences    % x y coordinates and ID (from topotoolbox) (nx3 matrix)
        Channelheads   % x y coordinates and ID (from topotoolbox)  (nx3 matrix)
        S              % STREAMobj topotoolbox
        Channels       % Stream divided into channels from one singular points (channelheads, confluences or outlet) to another
        hillslopes     % set of hillslopes (see class hillslope)
    end
    
    methods(Static)
        function obj=test
%             obj=watershed;
%             outlet_coord_BV6=[364778.7,6834822.7];
%             critic_drainage_area=40000;
%             obj=obj.analyze_hillslopes('MNT_PF_5m.tif',outlet_coord_BV6,critic_drainage_area);
%             folder_directory='C:\Users\Jean\Documents\Données\Synthetic\temp';
%             name_watershed='Pleine Fougeres'; 
%             obj.save_hillslopes(folder_directory,name_watershed);
            obj=watershed;
            outlet_coord_BV6=[364778.7,6834822.7];
            critic_drainage_area=40000;
            obj=obj.analyze_hillslopes('MNT_PF_5m.tif',outlet_coord_BV6,critic_drainage_area);
% % %                       obj=watershed;
% % %             outlet_coord_BV6=[351525,6790117];
% % %             critic_drainage_area=40000;
% % %             obj=obj.analyze_hillslopes('MNT_Vilaine_5m.tif',outlet_coord_BV6,critic_drainage_area);
% % %             folder_directory='C:\Users\Jean\Documents\Données\Synthetic\temp';
% % %             name_watershed='Upper Vilaine'; 
% % %             obj.save_hillslopes(folder_directory,name_watershed);
%             obj=watershed;
%             outlet_coord_Binic=[268863.41 ,6849415.51];
%             critic_drainage_area=40000;
%             obj=obj.analyze_hillslopes('Ic_BV2.tif',outlet_coord_Binic,critic_drainage_area);
%             folder_directory='C:\Users\Jean\Documents\Données\Synthetic\MorphologicInputs';
%             name_watershed='Ic at Binic'; 
%             obj.save_hillslopes(folder_directory,name_watershed);
%             name_watershed2 = name_watershed(~isspace(name_watershed));
%             save(strcat(folder_directory,'\watershed',name_watershed2,'.mat'),'obj');
%             obj=watershed;
%             outlet_coord_Kerb=[169565.30,6784657.50];
%             critic_drainage_area=10000;
%             obj=obj.analyze_hillslopes('MNT_Kerb_5m.tif',outlet_coord_Kerb,critic_drainage_area);
% %             folder_directory='C:\Users\Jean\Documents\Données\Synthetic\MorphologicInputs';
%             folder_directory='C:\Users\Jean\Documents\Données\Synthetic\temp';
%             name_watershed='Kerbernez extended'; 
%             obj.save_hillslopes(folder_directory,name_watershed);
%             name_watershed2 = name_watershed(~isspace(name_watershed));
%             save(strcat(folder_directory,'\watershed',name_watershed2,'.mat'),'obj');
            
%             obj=obj.load_DEM('Naizin_surrounding.tif');
%             outlet_coord_naizin=[266520.40,6779439.00];
%             obj=obj.get_watershed_DEM(outlet_coord_naizin);
% 
%             obj=obj.load_DEM('Blavet.tif');
%             outlet_coord_blavet=[255371.81,6774575.91];
%             obj=obj.get_watershed_DEM(outlet_coord_blavet);
%             obj=obj.load_DEM('Ic_BV2.tif');
%             outlet_coord_Binic=[268863.41 ,6849415.51];
%             obj=obj.get_watershed_DEM(outlet_coord_Binic);
%             obj=obj.load_DEM('MNT_Kerb_5m.tif');
%             outlet_coord_E30=[169091.7,6784372.8]; %E6=[118189.57,2347135.81];
%             obj=obj.get_watershed_DEM(outlet_coord_E30);
%             obj=obj.extract_stream_and_singular_points;
%             obj=obj.extract_channels_from_stream;
%             obj=obj.set_hillslopes_from_channels;
%             obj=obj.create_hillslopes1D;
%             folder_directory='C:\Users\Jean\Documents\Données\Synthetic\MorphologicInputs';
%             name_watershed='Pleine Fougeres'; 
%             name_watershed='Nazin extended'; 
%             folder_directory='C:\Users\Jean\Documents\Données\HillslopesData\Parameters';
%             folder_directory='C:\Users\Jean\Documents\Données\SyntheticData\MorphologicInputs';
%             obj.save_hillslopes(folder_directory,name_watershed);
        end
      
        function generate_hillslopes_from_dems
            dem_filenames={'MNT_PF_5m.tif','MNT_Kerb_5m.tif','Naizin_surrounding.tif','Blavet.tif','Ic_BV2.tif'};
            name_watersheds={'Pleine Fougeres','Kerbernez extended','Naizin extended','Blavet','Ic at Binic'};
            outlet_coords=[364778.7,6834822.7;169565.30,6784657.50;266520.40,6779439.00;255371.81,6774575.91;268863.41 ,6849415.51];
            critic_drainage_areas=[40000,10000,40000,20000,40000];
            folder_directory='C:\Users\Jean\Documents\Données\Synthetic\MorphologicInputs';
            
            for j=1:length(dem_filenames)
                obj=watershed;
                obj=obj.analyze_hillslopes(dem_filenames{j},outlet_coords(j,:),critic_drainage_areas(j));
                obj.save_hillslopes(folder_directory,name_watersheds{j});
                name_watershed =name_watersheds{j};
                name_watershed = name_watershed(~isspace(name_watershed));
                save(strcat(folder_directory,'\watershed',name_watershed,'.mat'),'obj');
            end
        end
    end
    
    methods(Access=public)
        % Constructor
        function obj=watershed
        end
        
        % Load DEM
        function obj=load_DEM(obj,filename)
            file_directory=which(filename);
            obj.DEM=GRIDobj(file_directory);
            % remove aberrant values into nan
            obj.DEM.Z(obj.DEM.Z<-10)=nan;
            obj.DEM=fillsinks(obj.DEM);
        end
        
        function obj=get_watershed_DEM(obj,manual_outlet_coord)
            if(nargin<2) manual_outlet_coord=[];    end
            FD=FLOWobj(obj.DEM);
            A=flowacc(FD);
            
            if(isempty(manual_outlet_coord))
                DB = drainagebasins(FD);
                DB = shufflelabel(DB);
                
                STATS = regionprops(DB.Z,'PixelIdxList','Area','Centroid');
                [~,IX] = max([STATS.Area]);
                obj.DEM.Z(DB.Z~=IX)=NaN;
            else
                W=A>10000;
                S=STREAMobj(FD,W);
                [xriv,yriv] = snap2stream(S,manual_outlet_coord(1),manual_outlet_coord(2));
                DB = drainagebasins(FD,xriv,yriv);
%                 Outlets_coord=streampoi(FD,W,'outlets','xy');
%                 Outlets_coord=Outlets_coord(find(min(((Outlets_coord(:,1)-manual_outlet_coord(1)).^2+(Outlets_coord(:,2)-manual_outlet_coord(2)).^2).^0.5)));
%                 DB = drainagebasins(FD,Outlets_coord(1),Outlets_coord(2));
                obj.DEM.Z(DB.Z==0)=NaN;
            end
        end
        

        function obj=extract_stream_and_singular_points(obj,critic_drainage_area)
            if(nargin<2) 
                % if not specified chose 40000 cells so for a 5m DEM it is 1km2
                critic_drainage_area=40000; 
            end
            FD=FLOWobj(obj.DEM);
            A=flowacc(FD);
            W = A>critic_drainage_area;
            obj.S = STREAMobj(FD,W);
            
            % Extract singular points
            obj.Channelheads=streampoi(FD,W,'channelheads','xy');
            obj.Outlet=streampoi(FD,W,'outlets','xy');
            obj.Confluences=streampoi(FD,W,'confluences','xy');
            % Get Topotoolbox ID
            obj.Channelheads=[obj.Channelheads,streampoi(FD,W,'channelheads','ix')];
            obj.Outlet=[obj.Outlet,streampoi(FD,W,'outlets','ix')];
            obj.Confluences=[obj.Confluences,streampoi(FD,W,'confluences','ix')];
        end
        
        function obj=extract_channels_from_stream(obj)
            % Store all singular points in a vector
            Singular_points=[obj.Channelheads(:,3);obj.Confluences(:,3)];
            Singular_points_coord=[obj.Channelheads(:,:);obj.Confluences(:,:)];
            
            % Find first point upstream of the confluences
            Confluences_coord_prec=[];
            Confluences_prec=[];
            for k=1:length(obj.Confluences(:,1))
                Id=find(obj.S.IXgrid==obj.Confluences(k,3));
                Precedent_nodes=obj.S.ix(obj.S.ixc==Id);
                Confluences_coord_prec=[Confluences_coord_prec; obj.S.x(Precedent_nodes),obj.S.y(Precedent_nodes)];
                Confluences_prec=[Confluences_prec; obj.S.IXgrid(Precedent_nodes)];
            end
            
            Singular_points_prec=[Confluences_prec;obj.Outlet(:,3)];
%             Singular_points_coord_prec=[Confluences_coord_prec;obj.Outlet(:,1:2)];
            
            obj.Channels=cell(1,length(Singular_points)+length(obj.Channelheads(:,1)));
            for i=1:length(Singular_points)
                obj.Channels{i}=Singular_points_coord(i,:);
                Test=[];
                A=Singular_points(i);
                while(isempty(Test))
                    Id=find(obj.S.IXgrid==A);
                    Next_node=obj.S.ixc(obj.S.ix==Id);
                    A=obj.S.IXgrid(Next_node);
                    Test=find(Singular_points_prec==A);
                    obj.Channels{i}=[obj.Channels{i};obj.S.x(Next_node),obj.S.y(Next_node),obj.S.IXgrid(Next_node)];
                end
            end
            
            for i=(length(Singular_points)+1):(length(Singular_points)+length(obj.Channelheads(:,1)))
                obj.Channels{i}=obj.Channelheads(i-length(Singular_points),1:3);
            end
        end
        
        function obj=set_hillslopes_from_channels(obj)
            FD=FLOWobj(obj.DEM);
            Distances_to_stream=flowdistance(FD,obj.S);
            compt=1;
            for i=1:length(obj.Channels)
                obj.hillslopes{i}=hillslope(compt);
                obj.hillslopes{i}=obj.hillslopes{i}.extract_hillslopes_from_channel(obj.DEM,obj.Channels{i},FD,Distances_to_stream);
                compt=compt+1;
                if(length(obj.Channels{i}(:,1))>1)
                    [obj.hillslopes{i}(1),obj.hillslopes{i}(2)]=obj.hillslopes{i}.extract_unique_hillslope_from_channel(obj.Channels{i});
                    compt=compt+1;
                end
            end
        end
        
        function obj=create_hillslopes1D(obj)
            for i=1:length(obj.hillslopes)
                if(length(obj.hillslopes{i})==1)
                    obj.hillslopes{i}=obj.hillslopes{i}.compute_hillslope1D; 
                else
                    obj.hillslopes{i}(1)=obj.hillslopes{i}(1).compute_hillslope1D;
                    obj.hillslopes{i}(2)=obj.hillslopes{i}(2).compute_hillslope1D;
                end
            end
        end
        
        function save_hillslopes(obj,folder_directory,name_watershed)
            for i=1:length(obj.hillslopes)
                if(length(obj.hillslopes{i})==1)
                    obj.hillslopes{i}.save_hillslope(folder_directory,name_watershed);
                else
                    obj.hillslopes{i}(1).save_hillslope(folder_directory,name_watershed);
                    obj.hillslopes{i}(2).save_hillslope(folder_directory,name_watershed);
                end
            end
        end
        
        function obj=analyze_hillslopes(obj,dem_filename,outlet_coord,critic_drainage_area)
            obj=obj.load_DEM(dem_filename);
            if(outlet_coord==-1)
                obj=obj.get_watershed_DEM;
            else
                obj=obj.get_watershed_DEM(outlet_coord);
            end
            if(critic_drainage_area==-1)
                obj=obj.extract_stream_and_singular_points;
            else
                obj=obj.extract_stream_and_singular_points(critic_drainage_area);
            end
            obj=obj.extract_channels_from_stream;
            obj=obj.set_hillslopes_from_channels;
            obj=obj.create_hillslopes1D;
        end
    end
end