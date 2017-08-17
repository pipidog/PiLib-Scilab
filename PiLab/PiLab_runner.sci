// Scilab command ====================================
function PiLab_runner(run)
    run.project_folder=PIL_dir_path(run.project_folder);
    cd(run.project_folder);
    PiLab(run.project_name,run.task_name,run.task_action)
    if run.force_exit=='on'
        exit
    end
endfunction
