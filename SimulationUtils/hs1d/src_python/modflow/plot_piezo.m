%%Initialization
close all;
clearvars;
% name = 'Laita';
% name = 'regular_watershed';
% name = 'Odet';
name = 'Kerbernez';

%% Import data from text file based on Watershed's name
path = pwd;
file_name = strcat(path,'/',name,'/',name,'_piezometry');
file_name_init = strcat(path, '/', name, '/', name, '_initial_state');
piezo = dlmread(file_name);
init = dlmread(file_name_init);
init = flipud(init);
piezo(piezo<-500) = NaN;
piezo(piezo>500) = NaN;
init(init<-500) = NaN;
init(init>500) = NaN;

%% Import X, Y and Z coordinates based on Watershed's name
path = strcat(pwd,'/',name,'/');
file_x = strcat(path,'coord_X_',name);
file_y = strcat(path,'coord_Y_',name);
file_z = strcat(path,name,'_elevation');

X_coord = dlmread(file_x);
X = unique(X_coord);

Y_coord = dlmread(file_y);
Y = unique(Y_coord);

Z_temp = dlmread(file_z);

Z = zeros(length(Y),length(X))-900;
for i = 1:length(Z_temp)
    a = find(X_coord(i) == X);
    b = find(Y_coord(i) == Y);
    Z(b,a) = Z_temp(i);
end

Z(Z<-500) = NaN;
piezo(piezo<Z-20)=NaN;

piezo = flipud(piezo);
Z = flipud(Z);

%% Plot the different surfaces

%Topographic level of the watershed
figure(1)
surface(X,Y,Z,'EdgeColor','none')
xlabel('X')
ylabel('y')
title(strcat(name,' elevation (m)'))
colorbar

figure(2)
surface(X,Y,piezo,'EdgeColor','none')
xlabel('X')
ylabel('Y')
title(strcat(name,' piezometric level (m)'))
colorbar

figure(3)
z_new = abs(max(max(Z-piezo)) - (Z-piezo));
surface(X,Y,z_new,'EdgeColor','none')
caxis([min(min(z_new)) max(max(z_new))])
xlabel('X')
ylabel('Y')
title(strcat(name,' - Piezo from cell bottom'))
colorbar

figure(4)
z_new = abs(max(max(Z-init)) - (Z-init));
surface(X,Y,z_new,'EdgeColor','none')
caxis([min(min(z_new)) max(max(z_new))])
xlabel('X')
ylabel('Y')
title(strcat(name,' - Initial State'))
colorbar

figure(5)
surface(X,Y,init,'EdgeColor','none')
xlabel('X')
ylabel('Y')
title(strcat(name,' initial piezometric level (m)'))
colorbar

figure(6)
z_new = piezo - init;
surface(X,Y,z_new,'EdgeColor','none')
caxis([min(min(z_new)) max(max(z_new))])
xlabel('X')
ylabel('Y')
title(strcat(name,' - Difference between Calc and init'))
colorbar


