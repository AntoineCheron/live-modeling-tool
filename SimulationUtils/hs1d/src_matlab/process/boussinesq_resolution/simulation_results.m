classdef simulation_results
    properties(Access=public)
        t
        x_Q
        x_S
        S
        Q
        QS
    end
    
    methods(Access=public)
        % constructor
        function obj=simulation_results(t,x_S,x_Q,S,Q,QS)
            obj.t=t;
            obj.x_S=x_S;
            obj.x_Q=x_Q;
            obj.S=S;
            obj.Q=Q;
            obj.QS=QS;
        end
        
        function plot_results(obj,S_max,Watershed_area,w,Flux_in,folder_results)
            % time parameters at which spatialized data are plotted
            tdays=time_unit.time_to_days(obj.t,'sec');
            tmax=tdays(end);
            t1=tmax/3; t2=3*tmax/4;
            t_initial=find((tdays-0).^2==min((tdays-0).^2));
            t_intermediar1=find((tdays-t1).^2==min((tdays-t1).^2));
            t_intermediar2=find((tdays-t2).^2==min((tdays-t2).^2));
            t_end=find((tdays-tmax).^2==min((tdays-tmax).^2));
            
            % compute relative storage matrix := S/Smax
            relative_storage=bsxfun (@rdivide, obj.S, S_max/100);
            
            % find which part of the hillslope is continuously saturated -> that is how we define the river location.
            Smin=min(relative_storage,[],2);
            pos=find(Smin==100);
            pos2=find(Smin<100);
            if(isempty(pos))
                river_pos=0;
                river_pos2=length(Smin)-1;
            else
                river_pos=pos(end);
                river_pos2=pos2(1)-1;
            end
            river_pos=1;
%             [r,c]=find(relative_storage==100);
%             pos=accumarray(c,r,[],@max);
%             idx_subsurface_flow= sub2ind(size(obj.Q),pos+1,(1:1:length(obj.Q))');
%             [r2,c2]=find(relative_storage<100);
%             idx_seepage_hillslope = sub2ind(size(obj.QS), r2, c2);
            % relative storage spatilized plot
            figure; hold on;
            plot(obj.x_S,relative_storage(:,t_initial(1)),'LineWidth',3,'LineStyle',':','Color',[0    0.4470    0.7410]);
            plot(obj.x_S,relative_storage(:,t_intermediar1(1)),'LineWidth',2,'LineStyle','-.','Color',[0.8500    0.3250    0.0980]);
            plot(obj.x_S,relative_storage(:,t_intermediar2(1)),'LineWidth',2,'LineStyle','--','Color',[0.4940    0.1840    0.5560]);
            plot(obj.x_S,relative_storage(:,t_end(end)),'LineWidth',3,'LineStyle','-','Color',[0.4660    0.6740    0.1880]);
            h=legend('Initial',[num2str(round(t1)),' days'],[num2str(round(t2)),' days'],[num2str(round(tmax)),' days']);
            xlabel('x [m]','FontSize',18,'FOntweight','bold');
            ylabel('S(x,t) / S_{c}(x,t) [%]','FontSize',18,'FOntweight','bold');
            title('Relative soil moisture storage','FontSize',20,'FOntweight','bold');
            set(h,'FontSize',17);
            box on;
            ylim([0 100]);
            set(gca,'linewidth',3,'fontsize',16,'fontweight','bold');
            if(nargin>5)
                set(gcf, 'Position', get(0, 'Screensize'));
                savefig([folder_results,'\storage_spatialized.fig']); 
                print([folder_results,'\storage_spatialized.png'],'-dpng');
                close all;
            end
            % subsurface flow spatialized plot
            figure; hold on;
            Q=obj.Q/Watershed_area*1e3*(24*3600);
            plot(obj.x_Q,-Q(:,t_initial(1)),'LineWidth',3,'LineStyle','-.');
            plot(obj.x_Q,-Q(:,t_intermediar1(1)),'LineWidth',2);
            plot(obj.x_Q,-Q(:,t_intermediar2(1)),'LineWidth',2);
            plot(obj.x_Q,-Q(:,t_end(end)),'LineWidth',3);
            h=legend('initial',[num2str(round(t1)),' days'],[num2str(round(t2)),' days'],[num2str(round(tmax)),' days']);
            xlabel('Distance [m]','FontSize',18,'FOntweight','bold');
            ylabel('Darcy Flux [mm/d]','FontSize',18,'FOntweight','bold');
            title('x-> Q(t,x)','FontSize',18,'FOntweight','bold');
            set(h,'FontSize',16);
            subsurface_values=-Q(river_pos+1:end,:);
            ymax=max(subsurface_values(:));
            ylim([0 ymax]);
            box on;
            set(gca,'linewidth',3,'fontsize',16,'fontweight','bold');
            if(nargin>5)
                set(gcf, 'Position', get(0, 'Screensize'));
                savefig([folder_results,'\subflow_spatialized.fig']);
                print([folder_results,'\subflow_spatialized.png'],'-dpng');
                close all;
            end
            % seepage flow spatialized plot
            figure; hold on;
            Seepage_local=diag(1./w)*obj.QS;               % convert it in m/s (obj.QS originally stored in m2/s (depending on the discretization chosen)
            Seepage_local=Seepage_local*1e3*(24*3600);     % convert it in mm/d
            plot(obj.x_S(2:end),Seepage_local(2:end,t_initial(1)),'LineWidth',3,'LineStyle',':','Color',[0    0.4470    0.7410]);
            plot(obj.x_S(2:end),Seepage_local(2:end,t_intermediar1(1)),'LineWidth',3,'LineStyle','-.','Color',[ 0.8500    0.3250    0.0980]);
            plot(obj.x_S(2:end),Seepage_local(2:end,t_intermediar2(1)),'LineWidth',3,'LineStyle','--','Color',[0.4940    0.1840    0.5560]);
            plot(obj.x_S(2:end),Seepage_local(2:end,t_end(end)),'LineWidth',3,'LineStyle','-','Color',[0.4660    0.6740    0.1880]);
            h=legend('Initial',[num2str(round(t1)),' days'],[num2str(round(t2)),' days'],[num2str(round(tmax)),' days']);
            xlabel('x [m]','FontSize',18,'FOntweight','bold');
            ylabel('q_{S}(x,t) [mm/d]','FontSize',18,'FOntweight','bold');
            title('Saturation excess overland flux','FontSize',20,'FOntweight','bold');
            set(h,'FontSize',17);
            box on;
            set(gca,'linewidth',3,'fontsize',16,'fontweight','bold');
            if(nargin>5)
                set(gcf, 'Position', get(0, 'Screensize'));
                savefig([folder_results,'\seepflow_spatialized.fig']);
                print([folder_results,'\seepflow_spatialized.png'],'-dpng');
                close all;
            end

            % some flux variables
            Seepage_spatialized=obj.compute_seepage_spatialized; % in m3/s (obj.QS originally stored in m2/s (depending on the discretization chosen)
            Outflow_river=sum(Seepage_spatialized,1); % total outflow flux (river+hillslope) in m3/s
            Outflow_river=Outflow_river/Watershed_area*1e3*24*3600; % in mm/d
            Seepage_flux_hillslope=sum(Seepage_spatialized(river_pos+1:end,:),1); % seepage flux on the hillslope only (m3/s)
            Seepage_flux_hillslope=Seepage_flux_hillslope/Watershed_area*1e3*24*3600; % in mm/d
            Subsurface_flux=-Q(river_pos+1,:); % subsurface flow in mm/d
            Flux_in=Flux_in/Watershed_area*1e3*24*3600;
            
%             Subsurface_flux=-Q(idx_subsurface_flow);
%             Seepage_flux_hillslope=Seepage_spatialized(idx_seepage_hillslope);
%             Seepage_flux_hillslope=accumarray(c2,Seepage_flux_hillslope,[],@nansum);
            % seepage & subsurface flow evolution
            figure; hold on;
            subplot(2,1,1);
            plot(tdays,Flux_in,'LineWidth',3); 
            ylabel('Flux normalized by Area [mm/d]','FontSize',18,'FOntweight','bold');
            h=legend('Infiltration');
            set(h,'FontSize',16);
            box on;
            set(gca,'linewidth',3,'fontsize',13,'fontweight','bold');
            subplot(2,1,2); hold on;
            plot(tdays,Outflow_river,'LineWidth',3);   % total outflow flux
            plot(tdays,Seepage_flux_hillslope,'LineWidth',3); % seepage flux on the hillslope
            plot(tdays,Subsurface_flux,'LineWidth',3);        % subsurface flux
            xlabel('Time [days]','FontSize',18,'FOntweight','bold');
            ylabel('Flux normalized by Area [mm/d]','FontSize',18,'FOntweight','bold');        
            h=legend('Outflow','Total seepage flux on the hillslope','Subsurface Flow');            
            set(h,'FontSize',16);
            box on;
            set(gca,'linewidth',3,'fontsize',13,'fontweight','bold');
            if(nargin>5)
                set(gcf, 'Position', get(0, 'Screensize'));
                savefig([folder_results,'\flow_vs_infiltration.fig']);
                print([folder_results,'\flow_vs_infiltration.png'],'-dpng');
                close all;
            end
            
            figure; hold on;
            plot(tdays,Outflow_river,'LineWidth',3);   % total outflow flux
            plot(tdays,Seepage_flux_hillslope,'LineWidth',3); % seepage flux on the hillslope
            plot(tdays,Subsurface_flux,'LineWidth',3);        % subsurface flux
            xlabel('Time [days]','FontSize',18,'FOntweight','bold');
            ylabel('Flux normalized by Area [mm/d]','FontSize',18,'FOntweight','bold');        
            h=legend('Outflow','Total seepage flux on the hillslope','Subsurface Flow');            
            set(h,'FontSize',16);
            box on;
            set(gca,'linewidth',3,'fontsize',13,'fontweight','bold');
            if(nargin>5)
                set(gcf, 'Position', get(0, 'Screensize'));
                savefig([folder_results,'\flow_evolution.fig']);
                print([folder_results,'\flow_evolution.png'],'-dpng');
                close all;
            end
        end
        
        function customed_plot_results(obj)
            figure;
            [hAx,hLine1,hLine2] = plotyy(tdays,Flux_in/Watershed_area*1e3*(24*3600),tdays,QSsumter);
            QSsum2=sum(QS,1);
            set(hAx,'NextPlot','add')
            plot(hAx(2),tdays,QSsum2)
%             plot(hAx(2),1:10,4:13)

            set(hLine1,'linewidth',2,'color','b','LineStyle','-')
            set(hLine2,'linewidth',2,'color','m','LineStyle','-')
            set(hAx(2),'xtick',[])
            set(get(hAx(1),'Ylabel'),'String','Rainfall [mm.h^{-1}]','fontweight','bold','fontsize',14)
            set(get(hAx(2),'Ylabel'),'String','Outflow [mm.d^{-1}]','fontweight','bold','fontsize',14)
            set(get(hAx(1),'Xlabel'),'String','Time [days]','fontweight','bold','fontsize',14)
            set(hAx,'FontSize',13,'fontweight','bold')
            set(hAx,{'ycolor'},{'b';'m'})
            
                figure; hold on;
            [hAx,hLine1,hLine2] = plotyy(tdays,Flux_in/Watershed_area*1e3*(24*3600),tdays,-Q(2,:));
            set(hLine1,'linewidth',2,'color','b','LineStyle','-')
            set(hLine2,'linewidth',2,'color','m','LineStyle','-')
            set(hAx(2),'xtick',[])
            set(get(hAx(1),'Ylabel'),'String','Rainfall [mm.h^{-1}]','fontweight','bold','fontsize',14)
            set(get(hAx(2),'Ylabel'),'String','Outflow [mm.d^{-1}]','fontweight','bold','fontsize',14)
            set(get(hAx(1),'Xlabel'),'String','Time [days]','fontweight','bold','fontsize',14)
            set(hAx,'FontSize',13,'fontweight','bold')
            set(hAx,{'ycolor'},{'b';'m'})
            
%             close all;
            figure;
            Q=obj.Q*1000;
            Q=-Q(3,:);
            QS=obj.compute_seepage_spatialized;
            QS=QS*1e3;
            QSsumter=Q+sum(QS(3,:),1);
            QSsum2=sum(QS,1);
            line(tdays,0.3*Flux_in*1e3*3600/Watershed_area,'Color',[0    0.4470    0.7410],'linewidth',3);
            legend('Infiltration')
            ax1 = gca;
            ax1_pos = ax1.Position; % position of first axes
            ax1.YColor = [0    0.4470    0.7410];
            ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
            ax2.YColor = [0.8500    0.3250    0.0980];
            line(tdays,QSsumter,'Parent',ax2,'Color',[0.8500    0.3250    0.0980],'linewidth',3)
            hold on
            line(tdays(1:3:end),QSsum2(1:3:end),'Parent',ax2,'Color',[0.9290    0.6940    0.1250],'LineStyle','none','Marker','x','linewidth',1)
%             legend('Subsurface Flow','Total Outflow (subsurface and seepage)')
            set(ax2,'xtick',[])
            set(get(ax1,'Ylabel'),'String','Infiltration [mm.h^{-1}]','fontweight','bold','fontsize',18)
            set(get(ax2,'Ylabel'),'String','Outflow [L.s^{-1}]','fontweight','bold','fontsize',18)
            set(get(ax1,'Xlabel'),'String','Time [days]','fontweight','bold','fontsize',18)
            set(ax1,'linewidth',3,'FontSize',16,'fontweight','bold')
            set(ax1,'ycolor',[0    0.4470    0.7410])
            set(ax2,'linewidth',3,'FontSize',16,'fontweight','bold')
            set(ax2,'ycolor',[0.8500    0.3250    0.0980])    
        end
        
        function save_results(obj,S_max,Watershed_area,foldername)
            filename_S=strcat(foldername,'\storage.code');
            fid = fopen(filename_S, 'w');
            fprintf(fid, 'First row: t, First column: x, S in the matrix (2:end,2:end) size(x,t) [m2] \n');
            fclose(fid);            
            S=nan(size(obj.S)+1);
            S(1,:)=[-1, obj.t];
            S(:,1)=[-1; obj.x_S];
            S(2:end,2:end)=obj.S;
            dlmwrite(filename_S,S, '-append', 'precision', '%E','delimiter','\t');
            
            filename_QS=strcat(foldername,'\seepage_flux.code');
            fid = fopen(filename_QS, 'w');
            fprintf(fid, 'First row: t, First column: x, QS in the matrix (2:end,2:end) size(x,t) [m2/s] \n');
            fclose(fid);         
            QS=nan(size(obj.QS)+1);
            QS(1,:)=[-1, obj.t];
            QS(:,1)=[-1; obj.x_S];
            QS(2:end,2:end)=obj.QS;
            dlmwrite(filename_QS,QS, '-append', 'precision', '%E','delimiter','\t');
            
            filename_Q=strcat(foldername,'\subsurface_flux.code');
            fid = fopen(filename_Q, 'w');
            fprintf(fid, 'First row: t, First column: x, Q in the matrix (2:end,2:end) size(x,t) [m3/s] \n');
            fclose(fid);
            Q=nan(size(obj.Q)+1);
            Q(1,:)=[-1, obj.t];
            Q(:,1)=[-1; obj.x_Q];
            Q(2:end,2:end)=obj.Q;
            dlmwrite(filename_Q,Q, '-append', 'precision', '%E','delimiter','\t');
            
            % compute relative storage matrix := S/Smax
            relative_storage=bsxfun (@rdivide, obj.S, S_max/100);
            filename_relstor=strcat(foldername,'\relative_storage.state');
            fid = fopen(filename_relstor, 'w');
            fprintf(fid, 'First row: t, First column: x, relative storage in the matrix (2:end,2:end) size(x,t) [percent] \n');
            fclose(fid);            
            relstor=nan(size(obj.S)+1);
            relstor(1,:)=[-1, obj.t];
            relstor(:,1)=[-1; obj.x_S];
            relstor(2:end,2:end)=relative_storage;
            dlmwrite(filename_relstor,relstor, '-append', 'precision', '%E','delimiter','\t');
            
            % find which part of the hillslope is continuously saturated -> that is how we define the river location.
            Smin=min(relative_storage,[],2);
            pos=find(Smin==100);
            if(isempty(pos))
                river_pos=0;
            else
                river_pos=pos(end);
            end
            
            filename_Q_river=strcat(foldername,'\river_flow_runoff.output');
            fid = fopen(filename_Q_river, 'w');
            fprintf(fid, 'First column: t, Q_river (total, at the channel) and runoff on the hillslope) in 2nd, 3rd and 4th column [m3/s]; columns 5th to 7th same metrics in [mm/d] \n');
            QS_out=obj.compute_seepage_spatialized;
            Q_river_total=sum(QS_out,1); Q_river_total=Q_river_total';
            if(river_pos==0)
                Q_river_channel=zeros(size(QS_out,2),1);
            else
                Q_river_channel=sum(QS_out(1:river_pos,:),1); Q_river_channel=Q_river_channel';
            end
            Q_seep_hillslope= Q_river_total-Q_river_channel;
            Q_river_total_mmd= Q_river_total/Watershed_area*1e3*24*3600;
            Q_river_channel_mmd= Q_river_channel/Watershed_area*1e3*24*3600;
            Q_seep_hillslope_mmd=Q_seep_hillslope/Watershed_area*1e3*24*3600;
            
            fclose(fid);
            Q_river=nan(length(obj.t),7);
            Q_river(:,1)=obj.t';
            Q_river(:,2)=Q_river_total;
            Q_river(:,3)=Q_river_channel;
            Q_river(:,4)=Q_seep_hillslope;
            Q_river(:,5)=Q_river_total_mmd;
            Q_river(:,6)=Q_river_channel_mmd;
            Q_river(:,7)=Q_seep_hillslope_mmd;
            dlmwrite(filename_Q_river,Q_river, '-append', 'precision', '%E','delimiter','\t');
            
            save(strcat(foldername,'\simulation_results.mat'),'obj');
        end
        
        function generate_storage_video(obj,S_max,w,Watershed_Area,Watershed_Area_spatialized,Flux_in,dS)
            
            fId = figure('Color','w');
            %             hold on
            % The final plot.
            %             plot();
            %% Set up the movie.
            writerObj = VideoWriter('C:\Users\Jean\Documents\temp\storage_video_N.mp4','MPEG-4'); % Name it.
            writerObj.FrameRate = 40; % How many frames per second.
            open(writerObj);
            Size_storage=length(obj.t);
            QS=obj.compute_seepage_spatialized;
            QS=diag(1./w)*obj.QS;
            QS=QS*1e3*(24*3600);
%             [Mass_balance,Flux_in,Darcy_Flux_out,Seepage_Flux_out,Storage_Var]=obj.compute_mass_changes(Watershed_Area,Flux_in,Watershed_Area_spatialized);
            [Mass_balance,Flux_in,Darcy_Flux_out,Seepage_Flux_out,Storage_Var]=obj.mass_changes_study(Watershed_Area_spatialized,Flux_in,dS);
            tdays=time_unit.time_to_days(obj.t,'sec');
            for i=1:Size_storage
                if mod(i,6)==0, %	At which frequence do you take frame for your video
                    figure(fId); % Makes sure you use your desired frame.
%                     subplot(1,3,1);
%                     set(gcf, 'Position', get(0, 'Screensize'));
%                     set(gcf, 'PaperPosition', [0 0 1280 1024]/200);
%                     set(gcf, 'Position', [353.8000  222.6000  966.4000  560.0000]);
%                     set(gcf, 'Position', [353.8000  222.6000  966.4000  620.0000]);
                    set(gcf, 'Position', [352  222  1024  640]);
%                     set(gcf, 'Position', [353.8000  222.6000  765  500.0000]);
%                     set(gcf, 'Position', [353.8000  222.6000  770.000  590.0000]);
                    
                    subplot(1,2,1)
                    hold off;
                    %             figure; hold on;
                    plot(tdays,Flux_in,'LineWidth',2);
                    hold on;
                    plot(tdays,Storage_Var,'LineWidth',2);
                    plot(tdays,-Darcy_Flux_out,'LineWidth',2);
                    plot(tdays,Seepage_Flux_out,'LineWidth',2);
                    plot(tdays,Mass_balance,'LineWidth',2);
                    h=legend('Infiltration','Storage variation','River discharge','Saturation excess overland flow','Total mass balance','Location','south');
                    xlabel('t [days]','FontSize',18,'FOntweight','bold');
                    ylabel('Flux [mm/d]','FontSize',18,'Fontweight','bold');
                    %             title('Mass balance','FontSize',18,'FOntweight','bold');
                    %             h2=title(['T= ',num2str(round(obj.t(i)/(24*3600))), ' days']);
                    
                    
                    box on;
                    set(gca,'linewidth',3,'fontsize',16,'fontweight','bold');
                    set(h,'FontSize',13,'fontweight','bold');
                    axis([0 35 -40 40]);
                    plot([obj.t(i)/(24*3600) obj.t(i)/(24*3600)],[-50 50],'r-','LineWidth',2);
                    hold off;
                    
                    subplot(1,2,2);
                    line(obj.x_S,obj.S(:,i)./S_max*100,'Color',[0    0.4470    0.7410],'linewidth',3);
                    ax1 = gca;
                    ax1_pos = ax1.Position; % position of first axes
                    ax1.YColor = [0    0.4470    0.7410];
                    ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
                    ax2.YColor = [0.8500    0.3250    0.0980];
                    line(obj.x_S(2:end),QS(2:end,i),'Parent',ax2,'Color',[0.8500    0.3250    0.0980],'linewidth',3,'LineStyle','-.')
                    set(ax2,'xtick',[])
                    set(get(ax1,'Ylabel'),'String','S(x,t) / S_{c}(x) [%]','fontweight','bold','fontsize',18)
                    set(get(ax2,'Ylabel'),'String','q_{S}(x,t) [mm.d^{-1}] (Dashed Line)','fontweight','bold','fontsize',18)
                    set(get(ax1,'Xlabel'),'String','x [m]','fontweight','bold','fontsize',18)
                    set(ax1,'linewidth',3,'FontSize',16,'fontweight','bold')
                    set(ax1,'ycolor',[0    0.4470    0.7410])
                    set(ax2,'linewidth',3,'FontSize',16,'fontweight','bold')
                    set(ax2,'ycolor',[0.8500    0.3250    0.0980])
                    xlim(ax1,[0 100]);
                    xlim(ax2,[0 100]);
                    ylim(ax1,[0 100]);
                    ylim(ax2,[0 80]);
                    
                    h2=suptitle(['Days simulated: ',num2str(round(obj.t(i)/(24*3600))), ' days']);
                    set(h2,'FontSize',16,'FontWeight','Bold');
                    
                    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                    writeVideo(writerObj, frame);
                end
                
            end
            hold off
            close(writerObj);
        end
        
        % m3/s -> integrate width function inside
        function QS_out_tot=compute_seepage_total(obj)
            QS_out=obj.compute_seepage_spatialized;
%             QS_out_tot=sum(QS_out);
            QS_out=QS_out(2:end,:);
            QS_out_tot=sum(QS_out);
        end
        
        % m3/s -> 
        function QS_out=compute_seepage_spatialized(obj)
            dx=obj.x_Q(2:end)-obj.x_Q(1:end-1);
            dx_QS=diag(dx);
            QS_out=dx_QS*obj.QS;
        end
        
        % m3/s -> how much is flowing out in the porous medium
        function Q=compute_darcy_flux(obj)
            Q=-obj.Q(2,:);
        end
        
        function Q_river=compute_river_flow(obj)
            Q_river=-obj.Q(2,:);
            Q_river=Q_river+obj.compute_seepage_total;
        end
        
        function [DPES_tot,DPES_spatialized]=compute_DPES(obj,Watershed_area_spatialized,Smax,Flux_in)
            relstor=bsxfun(@rdivide,obj.S,Smax);
            Bool_saturated=relstor>=1;
            Bool_saturated=Bool_saturated(2:end,:);
            Flux_in_spat=Watershed_area_spatialized(2:end)*Flux_in;
            DPES_spatialized=Flux_in_spat;
            DPES_spatialized(Bool_saturated==0)=0;
            DPES_tot=sum(DPES_spatialized);
        end
        
        function Q_diff_spatialized=compute_Q_difference_spatialized(obj)
            Q_diff_spatialized=obj.Q(1:end-1,:)-obj.Q(2:end,:);
        end
        
        function Storage_variation_spatialized=compute_storage_variation_spatialized(obj)
            dS=obj.S(:,3:end)-obj.S(:,1:end-2);
            dt=obj.t(3:end)-obj.t(1:end-2);
%                         Storage_variation_spatialized=dS*diag(1./dt);
            Storage_variation_spatialized=bsxfun (@rdivide, dS, dt);
        end
        
        function Storage_variation=compute_storage_variation(obj)
            Storage_variation_spatialized=obj.compute_storage_variation_spatialized;
            dx=obj.x_Q(2:end)-obj.x_Q(1:end-1);
            dx=diag(dx);
%             Storage_variation=sum(dx*Storage_variation_spatialized);
            Storage_variation_total=dx(2:end,2:end)*Storage_variation_spatialized(2:end,:);
            Storage_variation=sum(Storage_variation_total);
%             S=sum(dx*obj.S,1);
%             dS=S(3:end)-S(1:end-2);
%             dt=obj.t(3:end)-obj.t(1:end-2);
%             Storage_variation=dS./dt;
        end
        
        function [Mass_balance_spatialized,Flux_in_spatialized,Darcy_Flux_difference_spatialized,Seepage_Flux_spatialized,Storage_Var_spatialized]...
                =compute_mass_changes_spatialized(obj,Flux_in,Watershed_area_spatialized,dS)
            Flux_in_spatialized=Watershed_area_spatialized*Flux_in;
            Darcy_Flux_difference_spatialized=obj.compute_Q_difference_spatialized;
            Seepage_Flux_spatialized=obj.compute_seepage_spatialized;
            if(nargin<4)
                Storage_Var_spatialized=nan;
                Mass_balance_spatialized=nan;
            else
                Storage_Var_spatialized=dS;
                Mass_balance_spatialized=Flux_in_spatialized+Darcy_Flux_difference_spatialized-Seepage_Flux_spatialized-Storage_Var_spatialized;
            end
        end
        
        function [Mass_balance_tot_mmd,Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,Storage_Var_mmd]=mass_changes_study(obj,Watershed_Area_spatialized,Flux_in,dS)
            if(nargin<4) dS=nan; end
            [Mass_balance_spatialized,Flux_in_spatialized,Darcy_Flux_difference_spatialized,Seepage_Flux_spatialized,Storage_Var_spatialized]=obj.compute_mass_changes_spatialized(Flux_in,Watershed_Area_spatialized,dS);
            % from spatialized to the watershed scale mass balance
            Darcy_Flux_tot=obj.aggregate_spatialized_flux_at_watershed_scale(Darcy_Flux_difference_spatialized);
            Seepage_tot=obj.aggregate_spatialized_flux_at_watershed_scale(Seepage_Flux_spatialized);
            Flux_in_tot=obj.aggregate_spatialized_flux_at_watershed_scale(Flux_in_spatialized);
            if(nargin<4)
                Mass_balance_tot=nan;
                Storage_Var=nan;
            else
                Mass_balance_tot=obj.aggregate_spatialized_flux_at_watershed_scale(Mass_balance_spatialized);
                Storage_Var=obj.aggregate_spatialized_flux_at_watershed_scale(Storage_Var_spatialized);
            end
            
            
            % convert it in mmd
            Watershed_area_tot=sum(Watershed_Area_spatialized(2:end));
            Flux_in_tot_mmd=obj.convert_flux_balance_mmd(Flux_in_tot,Watershed_area_tot);
            Darcy_Flux_tot_mmd=obj.convert_flux_balance_mmd(Darcy_Flux_tot,Watershed_area_tot);
            Seepage_tot_mmd=obj.convert_flux_balance_mmd(Seepage_tot,Watershed_area_tot);
            if(nargin<4)
                Storage_Var_mmd=nan;
                Mass_balance_tot_mmd=nan;
            else
                Storage_Var_mmd=obj.convert_flux_balance_mmd(Storage_Var,Watershed_area_tot);
                Mass_balance_tot_mmd=obj.convert_flux_balance_mmd(Mass_balance_tot,Watershed_area_tot);
            end
            % plot it
%             obj.plot_mass_balance(Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,Storage_Var_mmd,Mass_balance_tot_mmd);
        end
        
        function plot_mass_balance(obj,Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,Storage_Var_mmd,Mass_balance_tot_mmd,folder_results)
            tdays=time_unit.time_to_days(obj.t,'sec');
            figure; hold on;
            plot(tdays,Flux_in_tot_mmd,'LineWidth',3);
            plot(tdays,Storage_Var_mmd,'LineWidth',3);
            plot(tdays,-Darcy_Flux_tot_mmd,'LineWidth',3);
            plot(tdays,Seepage_tot_mmd,'LineWidth',3);
            plot(tdays,Mass_balance_tot_mmd,'LineWidth',3);
            h=legend('Infiltration','Storage variation','River discharge','Saturation excess overland flow','Total mass balance','Location','south');
            xlabel('t [days]','FontSize',18,'FOntweight','bold');
            ylabel('Flux [mm/d]','FontSize',18,'FOntweight','bold');
            title('Mass balance components','FontSize',20,'FOntweight','bold');
            set(h,'FontSize',16);
            box on;
            axis([0 35 -40 40]);
            set(gca,'linewidth',3,'fontsize',17,'fontweight','bold');
%             axis([0 35 -50 50]);
            if(nargin>=7)
                set(gcf, 'Position', get(0, 'Screensize'));
                savefig([folder_results,'\mass_balance.fig']);
                print([folder_results,'\mass_balance.png'],'-dpng');
                close all;
            end
        end
        
        function flux_aggregated=aggregate_spatialized_flux_at_watershed_scale(obj,flux_spatialized)
            flux_spatialized=flux_spatialized(2:end,:);
            flux_aggregated=sum(flux_spatialized,1);
        end
        
        function Flux_tot_mmd=convert_flux_balance_mmd(obj,Flux_tot,Watershed_area_tot)
            Flux_tot_mmd=Flux_tot/Watershed_area_tot*1e3*24*3600;
        end
        function Flux_tot_mmd=convert_flux_balance_m3s(obj,Flux_tot,Watershed_area_tot)
            Flux_tot_mmd=Flux_tot*Watershed_area_tot/(1e3*24*3600);
        end
        
%         % all in mm/d
%         function [Mass_balance,Flux_in,Darcy_Flux_out,Seepage_Flux_out,Storage_Var]=compute_mass_changes(obj,Watershed_area,Flux_in,Watershed_area_spatialized)
%             Watershed_area_without_river=sum(Watershed_area_spatialized(2:end));
%             Flux_in=Flux_in(2:end-1);
%             Darcy_Flux_out=obj.compute_darcy_flux;
%             Darcy_Flux_out=Darcy_Flux_out(2:end-1);
%             Seepage_Flux_out=obj.compute_seepage_total;
%             Seepage_Flux_out=Seepage_Flux_out(2:end-1);
%             Storage_Var=obj.compute_storage_variation;
%             Flux_in=Flux_in/Watershed_area_without_river*1e3*(24*3600);
%             Darcy_Flux_out=Darcy_Flux_out/Watershed_area_without_river*1e3*(24*3600);
%             Seepage_Flux_out=Seepage_Flux_out/Watershed_area_without_river*1e3*(24*3600);
%             Storage_Var=Storage_Var/Watershed_area_without_river*1e3*(24*3600);
%             
%             Mass_balance=Flux_in-Darcy_Flux_out-Seepage_Flux_out-Storage_Var;
%         end
%         
%         function plot_mass_changes(obj,Watershed_Area,Flux_in,Watershed_area_spatialized,folder_results)
%             tdays=time_unit.time_to_days(obj.t,'sec');
%             [Mass_balance,Flux_in,Darcy_Flux_out,Seepage_Flux_out,Storage_Var]=obj.compute_mass_changes(Watershed_Area,Flux_in,Watershed_area_spatialized);
%             figure; hold on;
%             plot(tdays(2:end-1),Flux_in,'LineWidth',2);
%             plot(tdays(2:end-1),Mass_balance,'LineWidth',2);
%             plot(tdays(2:end-1),Darcy_Flux_out,'LineWidth',2);
%             plot(tdays(2:end-1),Storage_Var,'LineWidth',2);
%             plot(tdays(2:end-1),Seepage_Flux_out,'LineWidth',2);
%             h=legend('In Flux (Rain)','Total Mass Balance','Darcy flux exiting hillslope','Storage Variation','Seepage Flux');
%             xlabel('Time [days]','FontSize',16,'FOntweight','bold');
%             ylabel('Flux [mm/d]','FontSize',16,'FOntweight','bold');
%             title('Mass balance','FontSize',18,'FOntweight','bold');
%             set(h,'FontSize',13);
%             box on;
%             set(gca,'linewidth',3,'fontsize',16,'fontweight','bold');
%             axis([0 35 -50 50]);
%             if(nargin>=5)
%                 set(gcf, 'Position', get(0, 'Screensize'));
%                 savefig([folder_results,'\mass_balance.fig']);
%                 print([folder_results,'\mass_balance.png'],'-dpng');
%                 close all;
%             end
%         end
        

%         function time_travel_ref=sample_travel_time(obj)
        function [time_travel,time_ref,N_in_spatialized]=sample_travel_time(obj,N_in)
            % generate all possible combination of (x,t) to sample travel
            % times across all the space_time (x,t) space when it rains!
            time_sample=obj.t(N_in>0);
%             time_sample=time_sample(1:100:end);
            
            [A,B] = meshgrid(obj.x_S,time_sample);
            c=cat(2,A',B');
            x_t_sets=reshape(c,[],2);
            time_travel=nan(length(x_t_sets),1);
            time_ref=nan(length(x_t_sets),1);
            seepage_value=nan(length(x_t_sets),1);
            N_in_spatialized=interpn(obj.t,N_in,x_t_sets(:,2));
            
            v=obj.convert_flux_to_velocity;
            dx=obj.x_S(2:end)-obj.x_S(1:end-1);
            dt=diag(dx)*(1./v);
            % compute for each (x,t) pair the time it takes to reach the outlet
%             parfor i=1:length(x_t_sets)
            parfor i=1:length(x_t_sets)
                if(N_in_spatialized(i)>0)
%                 [time_travel_ref(i,1),time_travel_ref(i,2)]=obj.compute_travel_time(x_t_sets(i,2),x_t_sets(i,1),dt);
%                 [time_travel_ref(i,1),time_travel_ref(i,2)]=obj.compute_travel_time_double_hillslope(x_t_sets(i,2),x_t_sets(i,1),dt);
                 [time_travel(i),time_ref(i),seepage_value(i)]=obj.compute_travel_time(x_t_sets(i,2),x_t_sets(i,1),dt);
%                  [time_travel(i),time_ref(i)]=obj.compute_travel_time_double_hillslope(x_t_sets(i,2),x_t_sets(i,1),dt);
                end
            end
            
            [time_travel,time_ref]=obj.weight_by_seepage(time_travel,time_ref,seepage_value);
        end
        
        function [time_travel,t_ref,seepage_value]=compute_travel_time(obj,t_start,x_start,dt)
            Nx=find(obj.x_S==x_start)-1;
            Nt=find(obj.t==t_start);
            time_travel=0;
            t_ref=t_start;
            seepage_value=obj.QS(Nx+1,Nt);
            while(Nx~=0 || Nx==length(obj.x_S))
                dt_step=dt(Nx,Nt);
                deltat=abs(dt_step);
                time_travel=time_travel+deltat;
                t_ref=t_ref+deltat;
                if(dt_step>=0)
                    Nx=Nx-1;
                else
                    Nx=Nx+1;
                end
                [~,Nt]=min(abs(obj.t-t_ref));
                if(Nt==length(obj.t))
                    t_ref=nan;
                    time_travel=nan;
                end
            end
        end
        
        function [time_travel_real,time_ref_real]=weight_by_seepage(obj,time_travel,time_ref,seepage_value)
            time_travel_real=time_travel.*(seepage_value<=0);
            time_ref_real=time_ref-time_travel.*(seepage_value>0);
%             time_travel_real(seepage_value>0)=0;
%             time_ref_real(seepage_value>0)=time_ref(seepage_value>0)-time_travel(seepage_value>0);
        end
        
        function [time_travel,t_ref]=compute_travel_time_double_hillslope(obj,t_start,x_start,dt)
            dt2=-dt;
            Nx=find(obj.x_S==x_start);
            Nt=find(obj.t==t_start);
            time_travel=0;
            t_ref=t_start;
            if(Nx<length(obj.x_S)/2)
%                 Nx=Nx+1;
                while(Nx~=floor(length(obj.x_S)/2))
                    deltat=dt2(Nx,Nt);
                    time_travel=time_travel+deltat;
                    t_ref=t_ref+deltat;
                    Nx=Nx+1;
                    Nt=find(min(abs(obj.t-t_ref)));
                    if(Nt==length(obj.t))
                        t_ref=nan;
                        time_travel=nan;
                    end
                end
            elseif(Nx>length(obj.x_S)/2)
                Nx=Nx-1;
                while(Nx~=0)
                    deltat=dt(Nx,Nt);
                    time_travel=time_travel+deltat;
                    t_ref=t_ref+deltat;
                    Nx=Nx-1;
                    Nt=find(min(abs(obj.t-t_ref)));
                    if(Nt==length(obj.t))
                        t_ref=nan;
                        time_travel=nan;
                    end
                end
            end
        end
        
%         function [time_travel,t_ref]=compute_travel_time(obj,t_start,x_start,v)
%             time_travel=0;
%             t_ref=t_start;
%             % to be sure to begin at the center of one mesh
%             [~,x_ind]=min(abs(x_start-obj.x_S));
%             [~,t_ind]=min(abs(t_ref-obj.t));
%             
%             dx=obj.x_S(2:end)-obj.x_S(1:end-1);
%             iteration_nb=(x_ind-1);
%             for i=1:iteration_nb
%                 delta_t=dx(x_ind);
%                 AA=(v(x_ind,t_ind));
%                 time_travel=time_travel+delta_t;
%                 t_ref=t_ref+delta_t;
%                 [t_step,t_ind]=min(abs(t_ref-obj.t));
%                 if(t_ind==1)
%                     if(t_step>(obj.t(t_ind+1)-obj.t(t_ind)))
%                         time_travel=nan;
%                         t_ref=nan;
%                     end
%                 elseif(t_step>(obj.t(t_ind)-obj.t(t_ind-1)))
%                     time_travel=nan;
%                     t_ref=nan;
%                 end
%                 x_ind=x_ind-dx(x_ind);
%             end
%         end
%         
        function Q=interpolate_Q_meshcenter(obj)
            Q=nan(size(obj.S));
            for i=1:length(obj.Q(1,:))
                Q(:,i)=interpn(obj.x_Q,obj.Q(:,i),obj.x_S);
            end
        end
        
        function v=convert_flux_to_velocity(obj)
            S=nan(size(obj.Q));
            for i=1:length(obj.S(1,:))
                S(:,i)=interpn(obj.x_S,obj.S(:,i),obj.x_Q);
            end
            S=S(2:end-1,:);
            v=-obj.Q(2:end-1,:)./S;
%             Q=obj.interpolate_Q_meshcenter;
%             wS=diag(w)*obj.S;
%             v=-Q./wS;
        end
        
        function [tdec,QSdec,dQSdtdec,Qdec,dQdtdec]=get_decrease_parts2(obj,Flux_in,time_after_rain)
            t_sample=(obj.t(2:end)+obj.t(1:end-1))/2;
            Flux_in_sample=(Flux_in(2:end)+Flux_in(1:end-1))/2;
            QS_sample=(obj.QS(:,2:end)+obj.QS(:,1:end-1))/2;
            QS_sample=sum(QS_sample,1);
            QS_aggregated=sum(obj.QS,1);
            dt_sample=(obj.t(2:end)-obj.t(1:end-1));
            dQSdt_sample=(QS_aggregated(2:end)-QS_aggregated(1:end-1))./dt_sample;
            
            Q_sample=(-obj.Q(2,2:end)-obj.Q(2,1:end-1))/2;
            dQdt_sample=(-obj.Q(2,2:end)+obj.Q(2,1:end-1))./dt_sample;
            
%             Rain_event=t_sample(Flux_in_sample~=0);
            Rain_event=t_sample(Flux_in_sample>2e-08);
            Rain_event=[0,Rain_event,t_sample(end)];
            Time_between_2events=Rain_event(2:end)-Rain_event(1:end-1);
            max_Time=max(Time_between_2events);
            max_Time_steps=floor(max_Time/(min(dt_sample)));
            tdec=nan(length(Time_between_2events),max_Time_steps);
            QSdec=nan(length(Time_between_2events),max_Time_steps);
            dQSdtdec=nan(length(Time_between_2events),max_Time_steps);
            Qdec=nan(length(Time_between_2events),max_Time_steps);
            dQdtdec=nan(length(Time_between_2events),max_Time_steps);
            for i=1:length(Time_between_2events)
                if(Time_between_2events(i)>(time_after_rain*3600*24))
                    [~,tminstep]=min(abs(t_sample-(Rain_event(i)+360)));
                    [~,tmaxstep]=min(abs(t_sample-(Rain_event(i+1)-360)));
                    Length_curve=length(t_sample(tminstep:tmaxstep));
                    tdec(i,1:Length_curve)=t_sample(tminstep:tmaxstep);
                    QSdec(i,1:Length_curve)=QS_sample(tminstep:tmaxstep);
                    dQSdtdec(i,1:Length_curve)=dQSdt_sample(tminstep:tmaxstep);
                    Qdec(i,1:Length_curve)=Q_sample(tminstep:tmaxstep);
                    dQdtdec(i,1:Length_curve)=dQdt_sample(tminstep:tmaxstep);
                end
            end
            
            tdec=tdec(~all(isnan(tdec),2),:);
            QSdec=QSdec(~all(isnan(QSdec),2),:);
            dQSdtdec=dQSdtdec(~all(isnan(dQSdtdec),2),:);
            Qdec=Qdec(~all(isnan(Qdec),2),:);
            dQdtdec=dQdtdec(~all(isnan(dQdtdec),2),:);
            
        end
        
        
        function [coef,R2,tdec_interp,Qdec_interp,dQdec]=get_linear_regression_decrease_parts(obj,Flux_in,time_after_rain,t,Q)
            if(nargin<3) 
                time_after_rain=7200;
            end
            if(nargin<4)
                [tdec,Qdec]=obj.get_decrease_parts(Flux_in,time_after_rain);
            else
                [tdec,Qdec]=obj.get_decrease_parts(Flux_in,time_after_rain,t,Q);
            end
            [tdec_interp,dQdec,Qdec_interp]=simulation_results.get_derivative_decrease_parts(tdec,Qdec);
            coef=nan(length(dQdec),2);
            R2=nan(length(dQdec),1);
            for i=1:length(dQdec)
                [coef_temp,R2_temp]=simulation_results.linear_regression(Qdec_interp{i},dQdec{i});
                coef(i,1:2)=real(coef_temp)';
                R2(i,1)=real(R2_temp);
            end
        end
        
        function [tdec,Qdec]=get_decrease_parts(obj,Flux_in,time_after_rain,t,Q)
            if(nargin<3) 
                time_after_rain=7200;
            end
            if(nargin<4)
                t_study_rec=obj.t;
                Q_study_rec=-obj.Q(2,:);
            else
                t_study_rec=t;
                Q_study_rec=Q;
            end
            t_event=t_study_rec(Flux_in>0.5e-7);
            delta_t=t_event(2:end)-t_event(1:end-1);
            delta_t2=delta_t>(5*3600*24);
            t_event1=t_event([delta_t2,logical(0)]);
            t_event2=t_event([logical(0),delta_t2]);
            pos_1=nan(size(t_event1));
            pos_2=nan(size(t_event1));
            tdec=cell(size(t_event1));
            Qdec=cell(size(t_event1));
            for i=1:length(t_event1)
                [~,pos_1(i)]=min(abs(t_study_rec-(t_event1(i)+time_after_rain)));
                [~,pos_2(i)]=min(abs(t_study_rec-(t_event2(i)-3600*2)));
                tdec{i}=t_study_rec(pos_1(i):pos_2(i));
                Qdec{i}=Q_study_rec(pos_1(i):pos_2(i));
%                 smooth_function = fit(tdec{i}', Qdec{i}',  'smoothingspline', 'SmoothingParam', 1e-5);
%             	Qdec{i}=smooth_function(tdec{i}(1:10:end)');
%                 tdec{i}=tdec{i}(1:10:end)';
                Qdec{i} = smooth(tdec{i},Qdec{i},0.5,'rloess');
                Qdec{i} =Qdec{i}';
            end
            
        end
        
        % compute proportion of flow that have been Saturated Overland FLow
        function [part_overland_flow,part_subsurface_flow,mean_overland_flow,mean_subsurface_flow]=overland_flow_quantity(obj,Flux_in,Watershed_area_spatialized)
            [~,Flux_in_tot_mmd,Darcy_Flux_tot_mmd,Seepage_tot_mmd,~]=obj.mass_changes_study(Watershed_area_spatialized,Flux_in);
            Watershed_area_tot=sum(Watershed_area_spatialized(2:end));
            Flux_in_tot=obj.convert_flux_balance_m3s(Flux_in_tot_mmd,Watershed_area_tot);
            Darcy_Flux_tot=-obj.convert_flux_balance_m3s(Darcy_Flux_tot_mmd,Watershed_area_tot);
            Seepage_tot=obj.convert_flux_balance_m3s(Seepage_tot_mmd,Watershed_area_tot);
            
            dt=obj.t(2:end)-obj.t(1:end-1);
            Seepage_tot=(Seepage_tot(2:end)+Seepage_tot(1:end-1))/2;
            Darcy_Flux_tot=(Darcy_Flux_tot(2:end)+Darcy_Flux_tot(1:end-1))/2;
            Flux_in_tot=(Flux_in_tot(2:end)+Flux_in_tot(1:end-1))/2;
            Flux_in_total=sum(Flux_in_tot.*dt);
            qS_passed=sum(Seepage_tot.*dt);
            Q_passed=sum(Darcy_Flux_tot.*dt);
           
            part_overland_flow=qS_passed/(qS_passed+Q_passed);
            part_subsurface_flow=Q_passed/(qS_passed+Q_passed);
            part_to_storvar=(Flux_in_total-qS_passed-Q_passed)/Flux_in_total;
            qS_passed=obj.convert_flux_balance_mmd(qS_passed,Watershed_area_tot);
            Q_passed=obj.convert_flux_balance_mmd(Q_passed,Watershed_area_tot);
            mean_overland_flow=mean(qS_passed);
            mean_subsurface_flow=mean(Q_passed);
        end
        
        % compute relative time when the hillslope is at some locations fully saturated [%]
        % boolean 1 if sauration excess overland flow occured during the simulation
        function [proportion_era_seepage_front,boolean_seof]=importance_sat_zone(obj,S_max)
            boolean_seof=0;
            relstor=bsxfun(@rdivide,obj.S,S_max);
            relstor=relstor(2:end,:);
            saturated_area=relstor>=1;
            saturated_area_existence=sum(saturated_area);
            saturated_area_existence=(saturated_area_existence(1:end-1)+saturated_area_existence(2:end))/2;
            Time_step_seepage=saturated_area_existence>1;
            dt=obj.t(2:end)-obj.t(1:end-1);
            proportion_era_seepage_front=sum(Time_step_seepage.*dt)/(obj.t(end));
            if(proportion_era_seepage_front>0)
                boolean_seof=1;
            end
        end
        
        function [mean_humid_zone_length_adim,mean_humid_zone_area_adim,std_humid_zone_length_adim,std_humid_zone_area_adim]=compute_humid_zone_properties(obj,S_max,Watershed_area_spatialized)
            relstor=bsxfun(@rdivide,obj.S,S_max);
            relstor=relstor(2:end,:);
            saturated_area=relstor>=0.97;
            dx=obj.x_Q(2:end)-obj.x_Q(1:end-1);
            Watershed_area_spatialized=Watershed_area_spatialized(2:end);
            Area_sat_time=(Watershed_area_spatialized)'*saturated_area;
            Area_sat_time_adim=Area_sat_time/sum(Watershed_area_spatialized);
            Length_sat_time=((dx(2:end))')*saturated_area;
            Length_sat_time_adim=Length_sat_time/obj.x_Q(end);
            mean_humid_zone_length_adim=mean(Length_sat_time_adim);
            std_humid_zone_length_adim=std(Length_sat_time_adim);
            mean_humid_zone_area_adim=mean(Area_sat_time_adim);
            std_humid_zone_area_adim=std(Area_sat_time_adim);
        end
    end
    
    methods(Static)
        function [tdecdQ,dQdec,Qdec_interp]=get_derivative_decrease_parts(tdec,Qdec)
            tdecdQ=cell(1,length(tdec));
            dQdec=cell(1,length(tdec));
            Qdec_interp=cell(1,length(tdec));
            for i=1:length(tdec)
                tdecdQ{i}=(tdec{i}(2:end)+tdec{i}(1:end-1))/2/(3600*24);
                dQdec{i}=(Qdec{i}(2:end)-Qdec{i}(1:end-1))./((tdec{i}(2:end)-tdec{i}(1:end-1))/(3600*24));
                Qdec_interp{i}=(Qdec{i}(2:end)+Qdec{i}(1:end-1))/2;
            end
        end
        
        function [coef,R2]=linear_regression(Qdec,dQdec)
            X = [ones(length(log10(Qdec)),1) log10(Qdec)'];
            y=log10(-dQdec)';
            coef=X\y;
            yCalc = X*coef;
            R2 = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2);
        end
        
        function obj=load_sim_results_from_txt_files(sim_file_path)
            S_path=[sim_file_path,'\storage.code'];
            QS_path=[sim_file_path,'\seepage_flux.code'];
            Q_path=[sim_file_path,'\subsurface_flux.code'];
            [t_S,x_S,S]=simulation_results.read_storage_data(S_path);
            [t_Q,x_Q,Q]=simulation_results.read_storage_data(Q_path);
            [t_QS,x_QS,QS]=simulation_results.read_storage_data(QS_path);
            detect_diff_time=sum((t_S-t_Q).^2+(t_S-t_QS).^2);
            detect_diff_space=sum((x_S-x_QS).^2);
            if(detect_diff_time~=0)
                fprintf('WARNING: No coherent storages in simulation txtfiles for time sampling \n');
            end
            if(detect_diff_space~=0)
                fprintf('WARNING: No coherent storages in simulation txtfiles for space (x) sampling \n');
            end
            obj=simulation_results(t_S,x_S,x_Q,S,Q,QS);
        end
        
        function [t_S,x_S,S]=read_storage_data(S_path)
            fid_S=fopen(S_path);
            if(fid_S>0)
%                 textscan(fid_S, '%f', 'HeaderLines',1,'delimiter','\t' );
                 storage_data=dlmread(S_path,'\t' ,1,0);
            end
            t_S=storage_data(1,2:end);
            x_S=storage_data(2:end,1);
            S=storage_data(2:end,2:end);
            fclose(fid_S);
        end
        
        function [t_Q,x_Q,Q]=read_subsurface_data(Q_path)
            fid_Q=fopen(Q_path);
            if(fid_Q>0)
%                 textscan(fid_S, '%f', 'HeaderLines',1,'delimiter','\t' );
                 subsurface_data=dlmread(Q_path,'\t' ,1,0);
            end
            t_Q=subsurface_data(1,2:end);
            x_Q=subsurface_data(2:end,1);
            Q=subsurface_data(2:end,2:end);
            fclose(fid_Q);
        end
        
        function [t_QS,x_QS,QS]=read_seepage_data(QS_path)
            fid_QS=fopen(QS_path);
            if(fid_QS>0)
%                 textscan(fid_S, '%f', 'HeaderLines',1,'delimiter','\t' );
                 seepage_data=dlmread(QS_path,'\t' ,1,0);
            end
            t_QS=seepage_data(1,2:end);
            x_QS=seepage_data(2:end,1);
            QS=seepage_data(2:end,2:end);
            fclose(fid_QS);
        end
        
        function [t_S,x_S,S_rel]=read_relative_storage_data(S_path)
            fid_S=fopen(S_path);
            if(fid_S>0)
%                 textscan(fid_S, '%f', 'HeaderLines',1,'delimiter','\t' );
                 storage_data=dlmread(S_path,'\t' ,1,0);
            end
            t_S=storage_data(1,2:end);
            x_S=storage_data(2:end,1);
            S_rel=storage_data(2:end,2:end);
            fclose(fid_S);
        end
    end
end