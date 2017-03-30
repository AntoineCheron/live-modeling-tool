function folder_name=folder_of_file(file_name)
% Gets the containing folder
    folder_name=which(file_name); 
    a=strfind(folder_name,'\'); 
    folder_name=folder_name(1:a(end)); 
end