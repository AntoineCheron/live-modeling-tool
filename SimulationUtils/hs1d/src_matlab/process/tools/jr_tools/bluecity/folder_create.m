% Creates folder if does not exist
function folder_create(folder)
    if(exist(folder,'dir')~=7)
        success=mkdir(folder);
        if(success==0); fprintf('impossible to create directory %s\n',folder); end
    end
end
