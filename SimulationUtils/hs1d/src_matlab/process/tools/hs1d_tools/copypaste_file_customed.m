function destination_file=copypaste_file_customed(file_directory_source,folder_directory_destination)
    pos_file_name=strfind(file_directory_source,'\');
    file_name=file_directory_source(pos_file_name(end)+1:end);
    destination_file=fullfile(folder_directory_destination,file_name);
    copyfile(file_directory_source,destination_file);
end