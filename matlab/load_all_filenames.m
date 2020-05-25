function [ch1_filenames,folder_name] = load_all_filenames()


folder_name = uigetdir;
if folder_name== 0
    ch1_filenames=0;
    return
end
%read all files
tif_folder_name=strcat(folder_name,'\\*.tif');
tif_files=dir(tif_folder_name);
nTotFiles=length(tif_files);
marked_files=zeros(nTotFiles,1);
ch1_filenames={};

nFilePairs=0;

for nFile=1:length(tif_files)
    
    fn1=tif_files(nFile).name;
    fl1=length(fn1);
    fl1=fl1-4;
    ch1_filenames{nFile,1}=fn1(1:(fl1));
    
    %{    
    if(~marked_files(nFile))
        fn1=txt_files(nFile).name;
        fl1=length(fn1);
        ind1=fl1-ch1sufflen+1;
        fn1_first_part=fn1(1:(ind1-1));

        if(ind1>0)
            %if found ch1 suffix file
            if(strcmp(fn1(ind1:fl1),ch1suffix))
                %mark it as used
                marked_files(nFile)=1;
                %look for second channel
                for nFile2=1:length(txt_files)
                        if(~marked_files(nFile2))
                            fn2=txt_files(nFile2).name;
                            fl2=length(fn2);
                            ind2=fl2-ch2sufflen+1;

                            if(ind2>0)
                                %if found ch2 suffix file
                                fn2_first_part=fn2(1:(ind2-1));
                                if(strcmp(fn2(ind2:fl2),ch2suffix) && strcmp(fn1_first_part,fn2_first_part))
                                    %found second filename!
                                    marked_files(nFile2)=1;
                                    nFilePairs=nFilePairs+1;
                                    ch1_filenames{nFilePairs,1}=fn1;
                                    ch2_filenames{nFilePairs,1}=fn2;
                                    
                                end
                            end
                        end                        

                end
            end
        end
    end
%}
end


end