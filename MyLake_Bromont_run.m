% === function prep for MyLake run ===  %%%% 
%
% Module to input the test values of the parmaters
% Code checked by TSA, xx.03.2005
% Last modified by TSA, 15.08.2006 (Az replaced by In_Az 10.03.06; Possibility to have NaN in Global rad. series, 15.08.06)

function []= MyLake_Bromont_run(M_start,M_stop,lake_name,kz_N0,c_shelter, i_scv, i_sct, swa_b0, swa_b1, I_scDOC,k_Chl,k_POP,k_POC,k_DOP,k_DOC,k_pdesorb_a,k_pdesorb_b, enable_sediment,enable_river_inflow)
addpath(genpath("MyLake_v2_Vansjo"));
% Inputs:
%       M_start : Model start date [year, month, day]
%       M_stop : Model stop date [year, month, day]
%    
% Outputs:
%		

tic
disp('Started at:')
disp(datetime('now'));

m_start=[M_start, 1, 1]; %
m_stop=[M_stop, 12, 31];

save_initial_conditions = false; % save final concentrations as initial for the next run

name_of_scenario = sprintf('./IO/%s/input_%s.txt',lake_name,lake_name);
%name_of_init_file = sprintf('Inputs/%s/mylake_initial_concentrations.txt',lake_name);

% get basic value
[lake_params, sediment_params] = load_params();

% Niva results RMSD = 130 =======================================
file_name = sprintf('./Postproc_code/%s/%s_result_run.mat',lake_name,lake_name);
lake_params{47} = 58.3842e-003; % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
lake_params{49} = 128.2949e-003; % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
lake_params{50} = 1.4988e+000; % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
lake_params{53} = 1.6945e+000; % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
lake_params{56} = 208.3324e-003; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
lake_params{57} = 201.6135e-003; % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
lake_params{58} = 1.2687e+000; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
lake_params{59} = 1.6142e+000; % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
lake_params{46} = 31.3665e-003; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
lake_params{10} = 14.4699e-006; % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{54} = 30.5827e-006; % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{12} = 37.9560e-003; % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
lake_params{55} = 34.7141e-003; % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
sediment_params{52} = 21.5114e+000; % 65.1237e+000   %    accel
lake_params{24} = 373.1228e-003; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)


% Trials:
lake_params{24} = 1; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
lake_params{47} = 0.07; % 47     settling velocity for Chl1 a (m day-1)
lake_params{46} = 0.05; % 46  settling velocity for S (m day-1)
lake_params{56} = 0.07; % 56    Settling velocity for Chl2 a (m day-1)
sediment_params{52} = 100; % 65.1237e+000   %    accel
% =====================================================================================================================

% I_scO,k_BOD,I_scChl,
% Value giving during simulation

% ============ Water Column parameters ====================================
lake_params{4} = kz_N0;     %7.0e-05,     % 4     min. stability frequency (s-2)
if c_shelter == 'NaN'
    lake_params{5} = nan; %0.5,         % 5     wind shelter parameter (-)
else
    lake_params{5} = c_shelter; %0.5,         % 5     wind shelter parameter (-)
end
lake_params{16} = i_scv;    %1,           % 16    scaling factor for inflow volume (-)
lake_params{17} = i_sct;    %0,           % 17    adjusting delta for inflow temperature (-)
lake_params{39} = swa_b0;   %2.5,         % 39    non-PAR light attenuation coeff. (m-1)
lake_params{40} = swa_b1;   %1.05,        % 40    PAR light attenuation coeff. (m-1)
lake_params{23} = I_scDOC;   %1,           % 23    scaling factor for inflow concentration of DOC  (-)
% lake_params{25} = I_scO;    %1,             % 23    scaling factor for inflow concentration of DOC  (-)
% lake_params{22} = I_scChl;   %1           % 22    scaling factor for inflow concentration of Chl a (-)
% lake_params{22} = k_BOD;     %0.1         % 62    NOTE: NOT USED: Organic decomposition rate (1/d)
% =============== Sediment Parameters ====================================
%sediment_params{52} = 65.1237e+000;   %    accel
%lake_params{24} = 390.1162e-003;   % 24    scaling factor for inflow concentration of POP (-)
sediment_params{1} = k_Chl;                %        % 1
sediment_params{2} = k_POP;                %        % 1
sediment_params{3} = k_POC;                  %        % 0.01
sediment_params{4} = k_DOP;                 %        % 1
sediment_params{5} = k_DOC;                  %        % 1
sediment_params{23} = k_pdesorb_a; %100
sediment_params{24} = k_pdesorb_b; %100         % 

% for cores too (scaling unknown inputs):
% lake_params{18} = x(23);%    scaling factor for inflow concentration of C (-)
% lake_params{19} = x(24);%    scaling factor for inflow concentration of POC (-)
% lake_params{20} = x(25);%    scaling factor for inflow concentration of total P (-)
% lake_params{21} = x(26);%    scaling factor for inflow concentration of diss. organic P (-)
% lake_params{22} = x(27);%    scaling factor for inflow concentration of Chl a (-)
% lake_params{23} = x(28);%    scaling factor for inflow concentration of DOC  (-)
% lake_params{25} = x(29);%    Scaling factor for inflow concentration of O2 (-)
% lake_params{27} = x(30);%    Scaling factor for inflow concentration of NO3 (-)
% lake_params{34} = x(31);%    Scaling factor for inflow concentration of Fe3 (-)
% lake_params{35} = x(32);%    Scaling factor for inflow concentration of Al3 (-)
% lake_params{37} = x(33);%    Scaling factor for inflow concentration of CaCO3 (-)




% f = fopen('Inputs/vansjo_para.txt');
% garbage = fgetl(f); % file get line
% garbage = fgetl(f); % file get line
% data_lake = textscan(f, '%s%f%f%f%s', length(lake_params), 'Delimiter', '\t');
% fclose(f);
% 
lake_par_file = sprintf('./IO/%s/%s_para.txt',lake_name,lake_name);
if isfile(lake_par_file)
    
else
    export_params_lake(lake_params,'./IO/vansjo_para.txt', lake_par_file)
end

[MyLake_results, Sediment_results]  = fn_MyL_application_Bromont(m_start, m_stop, sediment_params, lake_params, name_of_scenario, save_initial_conditions, sprintf('./Postproc_code/%s',lake_name),sprintf('./IO/%s/mylake_initial_concentrations.txt',lake_name),lake_par_file, enable_sediment,enable_river_inflow); % runs the model and outputs obs and sim % runs the model and outputs obs and sim


disp('Saving results...')
save(file_name, 'MyLake_results', 'Sediment_results')
disp('Finished at:')
disp(datetime('now'));

toc
end
%
