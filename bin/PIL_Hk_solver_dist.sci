// **** Purpose ****
// This function can run k-point calculation distributedly. It split
// the k-point into several pieces and preared their input files for
// submission. You will need to assign how many tasks to submit each 
// time and the time interval between two batch submission.  
// So one can calculate a large k-mesh problem via a distirbutive way 
// to save time and memory.  
// **** Variables ****
// < Inputs >>
// [work_folder]: 1x1, string
// <= the folder where XXX_ham.sod file is stored. 
// [project_name]: 1x1, string
// <= project name of the ham file.
// [call_command]: 1 x 1, string 
// <= how to execute a file to scilab in your system. input and out
//    file names must be @input_name and @output_name. 
//    This code will automatically assign the names.
//    e.g. 'scilab-cli -nwni @input_name > @output_name &'  
// [k_point]: n x 3, real
// <= k-points to perform calculations
// [output_var]: 1x1, int
// <= output variable, 3: Ek, Vk, Hk; 2: Ek, Vk; 1: Ek
// [n_alloc]: 1 x 1, integer
// <= how many allocations of the k-points
// [n_thread]: 1x1, integer
// <= how many threads for each submission
// [thread_mem]: 1x1, real
// <= how much memory for each thread in unit of MB. Note that
//    thread_mem*n_thread should not more than the RAM safe limit 
//    of your machine.
// < Outputs >>
// [alloc_k_pt]: n_alloc*3, int
// => shows how the k-points are distributed. [tot_k,ini_k, final_k]
// 
// When running this function, it will also generate a folder called 
// "TaskName". Inside this folder, there will be several subfolder:
// [k_point]: k_point.sod, k_point_n.sod
// => k_point.sod stores all the k_points, k_point_n.sod store each piece
// of k_points.
// [k_input]: k_input_n.sce, n_input_n.log
// => input files and log files of each thread 
// [Ek], [Vk], [Hk]: eigenvalue, eigenvector, Hk 
// => Ek_n.sod: tot_state x tot_thread_n_k, eigenvelues
// => Vk_n.sod: tot_state x tot_state x tot_thread_n_k, eigenvectors
// => Hk_n.sod: tot_state x tot_state x tot_thread_n_k, Hamiltonian 
// There are also three files locate in the ./task_name:
// project_task.log shows detilas information of the calculaiton
// k_point_all.dat, k_point_all.sod shows all the k-points
// **** Version ****
// 03/28/2016 first built
// 05/18/2016 add job completion auto check. no time_int needed !
// **** Comment ****
// 1. To use this code, you must create a folder and put the ham file
//    in this folder as a working folder.
// 2. Once this code runs, this code will seperate the k-point into 
//    several segmenets and submit the job separately. When it is done,
//    it will generate a folder called 'k_dist'. Inside this folder,
//    there are several subfolder that contains the calculated results.
// 3. The output of k-point will automatically resize to n x 4 where the 
//    first coulmn shows it label. So k-points will not be messy in 
//    parallel approach. 
// 
function [alloc_k_pt]=PIL_Hk_solver_dist(work_folder,project_name,..
    task_name,call_command,k_point_all,output_var,n_alloc,n_thread,thread_mem)
    // thread_mem check
    if exists('thread_mem')==0 then
        thread_mem=1e7;
    end
    if thread_mem>2040.0 then
        thread_mem=2040.0;
        disp('Warning: thread_mem is too high, reset to 2040.0 !');
    end

    // correct work_folder name
    work_folder=PIL_dir_path(work_folder);
    task_folder=work_folder+task_name+'/';
    mkdir(task_folder);

    // check call command
    if (grep(call_command,'@input_name')+grep(call_command,'@output_name'))~=2 then
        disp('Error: PIL_Hk_solver_dist, @input_name and @output_name'..
        +' must appear in call_command');
        abort;
    end

    // reconstruct call command
    // full command is achieved by setting call_command(2) and call_command(4)
    tmp=strsplit(call_command);
    fn_index=find(tmp=='@');
    call_command=[strcat(tmp(1:fn_index(1)-1)),'',..
    strcat(tmp(fn_index(1)+11:fn_index(2)-1)),'',..
    strcat(tmp(fn_index(2)+12:$))]

    // generate parallel input files    
    fid(1)=mopen(task_folder+'runtime.log','w');
    mfprintf(fid(1),'    PIL_Hk_solver_dist: Creating folders\n\n');
    mfprintf(fid(1),'    PIL_Hk_solver_dist: Requesting memory \n')
    mfprintf(fid(1),'      Memory per allocation=%f MB\n\n',thread_mem)

    mkdir(task_folder+'k_point');
    mdelete(task_folder+'k_point/*.sod');

    mkdir(task_folder+'k_input');
    mdelete(task_folder+'k_input/*.sce');
    mdelete(task_folder+'k_input/*.log');

    mkdir(task_folder+'Ek');
    mdelete(task_folder+'Ek/*.sod');

    if output_var>=2  then
        mkdir(task_folder+'Vk');
        mdelete(task_folder+'Vk/*.sod');
    end

    if output_var>=3 then
        mkdir(task_folder+'Hk');
        mdelete(task_folder+'Hk/*.sod');
    end

    // output all k-point
    fid(2)=mopen(task_folder+'k_point_all.dat','w');
    PIL_print_mat('k_point_all, @f:f, all k-point',..
    k_point_all,'r',fid(2),'off');
    mclose(fid(2));

    // generate k distrubution variable
    tot_k=length(k_point_all(:,1));
    n_k=ceil(tot_k/n_alloc);  // num of k-pt except the last task
    alloc_k_pt=zeros(n_alloc,3);
    alloc_k_pt(:,1)=n_k*ones(n_alloc,1);
    alloc_k_pt($,1)=tot_k-n_k*(n_alloc-1);
    count=0;
    for n=1:n_alloc
        alloc_k_pt(n,2)=count+1;
        alloc_k_pt(n,3)=count+alloc_k_pt(n,1);
        count=count+alloc_k_pt(n,1);
    end

    if find(alloc_k_pt(:,1)<=0)~=[] then
        disp('Error: PIL_Hk_solver_dist, negative allocation!'..
        +' Use less allocatlions or make tot_k/n_alloc=integer.');
        abort
    end

    // print initial message
    mfprintf(fid(1),..
    '    PIL_HK_solver_dist: Starting distributed k-point calculations\n');
    mfprintf(fid(1),'      running machine= %s\n',unix_g('hostname'));
    mfprintf(fid(1),'      total k-points= %d\n',tot_k);
    mfprintf(fid(1),'      total allocations= %d\n',n_alloc);
    mfprintf(fid(1),'      total threads= %d\n',n_thread);
    mfprintf(fid(1),'      total batch submission= %d\n\n',ceil(n_alloc/n_thread));
    mfprintf(fid(1),'      allocation %4d, number of k-points= %4d,'+..
    ' range=%7d ~ %7d\n',[1:n_alloc]',alloc_k_pt);
    mfprintf(fid(1),'\n');

    // preparing input files
    mfprintf(fid(1),'    PIL_Hk_solver_dist: preparing input files\n\n')
    save(task_folder+'k_point_all.sod','k_point_all');
    n_k=ceil(tot_k/n_alloc);
    for n=1:n_alloc
        if n~=n_alloc
            k_point=k_point_all((n-1)*n_k+1:n*n_k,:);
        else
            k_point=k_point_all((n-1)*n_k+1:$,:);
        end
        save(task_folder+'k_point/k_point_'..
        +string(n)+'.sod','k_point')

        fid(2)=mopen(task_folder+'k_input/k_input_'..
        +string(n)+'.sce','w');
        mputl('clear; clc; stacksize('+string(fix(thread_mem/76*1e7))+..
        ') ;exec(PiLib);tic()',fid(2));
        
        mputl('printf(''\n\n'')',fid(2));
        mputl('printf(''    current time= %d/%d/%d %d:%d:%d\n'',clock())',fid(2));
        mputl('printf(''    Starting Hk solver on '+string(n)..
        +'-th allocation'')',fid(2))
        mputl('printf(''\n\n'')',fid(2))
        mputl('',fid(2));
        mputl('chdir('''+task_folder+''')',fid(2))
        mputl('load(''../'+project_name+'_ham.sod'')',fid(2));
        mputl('load(''./k_point/k_point_'+string(n)+'.sod'')',fid(2));
        select output_var
        case 3
            mputl('[Ek,Vk,Hk]=PIL_Hk_solver(k_point'..
            +',ham.lat_vec,ham.Hr_mat,ham.uc_index);',fid(2));
            mputl('save(''./Ek/Ek_'+string(n)+'.sod'',''Ek'')',fid(2));
            mputl('save(''./Vk/Vk_'+string(n)+'.sod'',''Vk'')',fid(2));
            mputl('save(''./Hk/Hk_'+string(n)+'.sod'',''Hk'')',fid(2));
        case 2
            mputl('[Ek,Vk]=PIL_Hk_solver(k_point'..
            +',ham.lat_vec,ham.Hr_mat,ham.uc_index);',fid);
            mputl('save(''./Ek/Ek_'+string(n)+'.sod'',''Ek'')',fid(2));
            mputl('save(''./Vk/Vk_'+string(n)+'.sod'',''Vk'')',fid(2));
        case 1
            mputl('[Ek]=PIL_Hk_solver(k_point'..
            +',ham.lat_vec,ham.Hr_mat,ham.uc_index);',fid(2));
            mputl('save(''./Ek/Ek_'+string(n)+'.sod'',''Ek'')',fid(2));
        end
        mputl('',fid(2));
        mputl('printf(''\n'')',fid(2));
        mputl('printf(''    time elapse=%f\n'',toc())',fid(2));
        mputl('printf(''    finish time= %d/%d/%d %d:%d:%d\n'',clock())',fid(2));
        mputl('printf(''\n\n'')',fid(2));
        mputl('exit()',fid(2));
        mclose(fid(2));
    end

    // submit parallel input files
    mfprintf(fid(1),'    PIL_Hk_solver_dist: batch job submission\n')
    mfprintf(fid(1),'\n')
    count=0;
    for n=1:n_alloc
        // input file
        call_command(2)=task_folder+'k_input/k_input_'+string(n)+'.sce';
        // output file
        call_command(4)=task_folder+'k_input/k_input_'+string(n)+'.log';
        // count batch submit
        if pmodulo(n,n_thread)==1 then
            count=count+1;
            mfprintf(fid(1),'      ## Batch Submit-%d\n\n',count);
        end

        // submit jobs 
        unix(strcat(call_command));
        mfprintf(fid(1),'      allocation %4d submitted\n',n)
        // show information
        if pmodulo(n,n_thread)==0 & n~=n_alloc
            c_tmp=clock();
            mfprintf(fid(1),'\n');
            mfprintf(fid(1),'      waiting for job completion\n');
            mfprintf(fid(1),'      current time= %4d/%02d/%02d %02d:%02d:%02d\n',c_tmp);
            mfprintf(fid(1),'\n');

            // check last run
            Ek_check=0;
            while Ek_check==0
                Ek_dir=dir(task_folder+'Ek');
                Ek_check=1;
                for m=1:n_thread
                    if grep(Ek_dir(2),'Ek_'+string(n-n_thread+m))==[] then
                        sleep(20*1e+3);
                        Ek_check=0;
                        break;
                    end
                end
            end
            mfprintf(fid(1),'      all jobs completed\n');
            mfprintf(fid(1),'      time elapse= %f seconds\n',etime(clock(),c_tmp));           
            mfprintf(fid(1),'\n');
        end 
    end
    mfprintf(fid(1),'    PIL_Hk_solver_dist: job submit completed\n\n');
    mclose(fid(1));
endfunction
