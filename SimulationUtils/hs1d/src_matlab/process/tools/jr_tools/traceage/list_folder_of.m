function list=list_folder_of(folder)
    % Lists all the subfolder of "folder" except . and ..: alphabetic order % order by creation time 
    temp=dir(folder); 
    [~,idx] = sort([temp.datenum]);
    temp=temp(idx);
    comp=1; 
    for i=1:length(temp)
        a=temp(i).name; 
        if(a(1)~='.')
            list{comp}=temp(i).name;
            comp=comp+1; 
        end
    end
    
end