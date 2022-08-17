function export_params_lake(lake_params,parfileplace, parfile)
    fileID = fopen(parfile,'w');
    fwrite(fileID,magic(5));
    fclose(fileID);
    basicfile = parfileplace;    
    if exist(basicfile, 'file') == 2
                % Read txt into cell A
        fid = fopen(basicfile,'r');
        i = 1;
        tline = fgetl(fid);
        A{i} = tline;
        while ischar(tline)
            i = i+1;
      
            if( i>2 )
                j = i-2;
                tline = fgets(fid);
                if length(tline) == 1
                    disp('end parfile');
                else
                    sc = split(tline);
                    cell = lake_params(j,1);
                    t0 = num2str(cell{1,1});
                    t1=sc{3};
                    t2=num2str(sc{4});
                    t3=num2str(sc{5});
                    t4=sc{6};
                    tline = sprintf('%s\t%s\t%s\t%s\t%s',sc{1},num2str(cell{1,1}),num2str(sc{3}),num2str(sc{4}),sc{5},sc{6});
                end
            else
                j = i-2;
                tline = fgets(fid);
                if j >= 61
                    tline = fgets(fid);
                else    
                    try
                    sc = split(tline);
                    catch
                        allo=1;
                    end
                    if j == 0
                      tline = sprintf('            %s\t%s\t%s\t%s\t%s',sc{2},sc{3},sc{4},sc{5},sc{6});

                    else
                        tline = sprintf('            %s\t%s\t%s\t%s\t%s',sc{2},num2str(sc{3}),num2str(sc{4}),num2str(sc{5}),sc{6});
                    end
                end
            end
            A{i} = tline;
        end
        fclose(fid);
        
        fid = fopen(parfile, 'w');
        for i = 1:numel(A)
           
                if A{i+1} == -1
                    fprintf(fid,'%s', A{i});
                    break
                else
                    fprintf(fid,'%s\n', A{i});
                end
            
        end 
    else
                % Read txt into cell A
        fid = fopen(parfile,'r');
        i = 1;
        tline = fgetl(fid);
        A{i} = tline;
        while ischar(tline)
            i = i+1;
            if i == 25 || i== 27 || i == 49|| i == 50
                j = i-2;
                tline = fgets(fid);
                sc = split(tline);
                cell = lake_params(j,1);
                tline = [sc{1},num2str(cell{1}),num2str(sc{3}),num2str(sc{4}),sc{5}];
            else
            tline = fgetl(fid);
            end
            A{i} = tline;
        end
        fclose(fid);
        fid = fopen(parfile, 'w');
        for i = 1:numel(A)
            if A{i+1} == -1
                fprintf(fid,'%s', A{i});
                break
            else
                fprintf(fid,'%s\n', A{i});
            end
        end 
         
    end
end
