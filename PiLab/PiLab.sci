// **** Purpose ****
// This is the PiLab caller. It provides a simple command to execute PiLab 
// functions without memorizing the all functions names of PiLab.
// **** Variables ****
// [project_name]: 1x1, string
// <= the project name
// [task_name]: 1x1, string
// <= the task_name to do 
// [action]: 1x1, string, default:'run'
// <= 'run', 'template', 'reload'
// **** Version ****
// 05/01/2014 first version
// 05/12/2015 add new keyword: action
// 05/22/2015 move input file loader in this file. 
// **** Comment ****
// action:
// 'run' : execute the task_name
// 'template' : generate the task_name
// 'reload' : re-generate .sod file based on the task_name text-file
function PiLab(project_name,task_name,action)
    [lhs,rhs]=argn();
    select rhs
    case 0
        project_name=input('project name = ','string');
        task_name=input('project task_name = ','string');
        action=input('action = ','string') 
    case 1
        disp('Error: PiLab, input arguments not enough!');
        abort
    case 2
        action='run'
    end
    select action
    case 'run'
        execstr('PiLab_'+task_name+'('''+project_name+''')')
    case 'template'
        PiLab_template(project_name,task_name);
    case 'reload'
        PiLab_loader(project_name,task_name,'all','keep');
    else
        disp('Error: PiLab, no corresponding action!');
        abort
    end
endfunction


