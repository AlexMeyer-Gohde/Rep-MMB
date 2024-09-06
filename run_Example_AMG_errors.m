% This file generate the results for the paper from the replication
% database from macromodeldatabase.com
% You can set run_worked_model=1 to run again with only worked model 
% after the first run through all model replications.


clear;


fprintf('-- %s.m started. --\n',mfilename);


% Inputs -------------------------------

folder_routine = 'Example_AMG_errors';
folder_routine_exec = 'Example_AMG_errors';
replication_directory_name = 'Rep-MMB';
overview_table_input_file = 'overview_out.xlsx';
final_file_all='Result_allmodels.xlsx';
final_file_worked='Result_worked_models';
matfile_all='AMG_Results.mat';
matfile_worked='AMG_Results_worked.mat';
% Inputs Ende --------------------------

% % create a test environment for the research routine
% prepare_folder(folder_routine_exec);
% % copy routine into test folder
% copyfile(folder_routine,folder_routine_exec,'f');
% % copy replication directory into test folder
% %copyfile(replication_directory_name, [folder_routine_exec '\' replication_directory_name], 'f');
% % copy overview table into test folder
% %copyfile(overview_table_input_file, [folder_routine_exec '\' overview_table_input_file]);
% toc
% fprintf('\n-- %s.m completed. --\n',mfilename);

cd(folder_routine)
mmb_loop;
%%%%%Choose to run again with only worked models
run_worked_model=1; % 1 if you want to run only model with Newton, 0 if you want to run the whole replication database 

if run_worked_model==1
mmb_loop_worked_models;
cd ..
folder_routine = 'Example_AMG_errors';
folder_routine_exec = 'Example_AMG_errors';
replication_directory_name = 'Rep-MMB';
overview_table_input_file = 'overview_out.xlsx';
final_file_all='Result_allmodels.xlsx';
final_file_worked='Result_worked_models';
matfile_all='AMG_Results.mat';
matfile_worked='AMG_Results_worked.mat';
%copyfile( [folder_routine_exec '\Results\' final_file_all]);   
%copyfile( [folder_routine_exec '\Results\' final_file_worked]);   
%copyfile( [folder_routine_exec '\Results\' matfile_worked]); 
copyfile( [folder_routine_exec '\Results\' matfile_all]); 
else
cd ..
folder_routine = 'Example_AMG_errors';
folder_routine_exec = 'Example_AMG_errors';
replication_directory_name = 'Rep-MMB';
overview_table_input_file = 'overview_out.xlsx';
final_file_all='Result_allmodels.xlsx';
final_file_worked='Result_worked_models';
matfile_all='AMG_Results.mat';
matfile_worked='AMG_Results_worked.mat';
%copyfile( [folder_routine_exec '\Results\' final_file_all]);   
copyfile( [folder_routine_exec '\Results\' matfile_all]);  
end


% function prepare_folder(folder)
% %PREPARE_FOLDER makes sure a folder of given path and name exists and is
% %empty.
% %   folder is a character array containing (relative) path and name of the
% %   folder to be created.
%     if isfolder(folder)
%         rmdir(folder,'s');
%     end
%     mkdir(folder);
% end