function [MyLake_results, Sediment_results] = fn_MyL_application_Bromont(m_start,m_stop, K_sediments, K_lake, name_of_scenario, is_save_results,~,initfile, param_file,  enable_sediment,enable_river_inflow)
global sed_par_file lake_par_file Eevapor
% This is the main MyLake application configuration file. INCA is a switch
% It is made to run a after the parameter are set by Set_Prior

Eevapor=0;
% disp('init ...');

calibration_k_values = [(1:length(K_sediments))',cell2mat(K_sediments(:,1)) ]; % writing sediments parameters file

%% generates unique files

sed_par_file = tempname;
lake_par_file = tempname;

dlmwrite(sed_par_file, calibration_k_values,'delimiter','\t');

%% writing lake parameter file
f = fopen(param_file);
garbage = fgetl(f); % file get line
garbage = fgetl(f); % file get line
data_lake = textscan(f, '%s%f%f%f%s', length(K_lake), 'Delimiter', '\t');
fclose(f); % the parameter line (xx,1) + 2 lines gives the location of the paramter in the input txt file.
% array position + 2 = input file line

% for i=1:length(K_lake)
%     data_lake{1, 2}(i,1) = K_lake{i}; % I_scDOC
% end


fid=fopen(lake_par_file,'wt');
fprintf(fid,'\n\n');
new_param = K_lake(:,1);
dlmwrite(lake_par_file, [[1:length(K_lake)]',cell2mat(new_param),data_lake{3},data_lake{4},(1:length(K_lake))'],'delimiter','\t','-append'); % 1:length(K_lake) is the length of the parameter file.
fclose(fid);
%% Specific MyLake application

% warning('off', 'all')

parafile=lake_par_file;
inputfile=name_of_scenario;

[MyLake_results, Sediment_results] ...
    = solvemodel_v2_Bromont(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake',enable_sediment,enable_river_inflow);

%MyLake_results = MyLake_results;
%Sediment_results = sediment_results;

% f1_name = (strcat(output_file, '/Tzt.csv')); % b = binary mode, z = archived file
% %f1_name = '\Tzt.csv';
% dlmwrite(f1_name, Tzt(:, 366:end)', 'delimiter', ',', 'precision', '%.3f'); % depend on spin up year
% O2zt = O2zt* 0.001;
% f5_name = (strcat(output_file, '/O2zt.csv'));
% dlmwrite(f5_name, O2zt(:, 366:end)', 'delimiter', ',', 'precision', '%.3f');
% 
% f5_name = (strcat(output_file, '/Chlzt.csv'));
% dlmwrite(f5_name,Chlzt(:, 366:end)', 'delimiter', ',', 'precision', '%.3f');

if is_save_results
    % disp('Saving sediment and water-column profiles for basin 1: Storefjorden');
    if enable_sediment == 1
        sediment_save_init_conc(Sediment_results, 1)
    end
    MyLake_save_result_for_init_conc(MyLake_results, 1)
else
    disp('Skipping saving the results and initial concentrations');
end


%% cleaning
fclose('all');
% delete (sed_par_file)
% delete (lake_par_file)

end
