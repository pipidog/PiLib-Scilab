// **** PURPOSE **** 
// This code generates the Fock state for a given occupation number
// **** VARIABLES ****
// [occup_num]: Nx2, integer
// <= tells the code the total state and occupation number of each site
// [index_form]: 1x1, string, 'total' or 'site'
// <= tells the code how will you label your Fock index 
//    (by total state or by each site) 
// [Fock_index]: MxP, integer, M: the total Fock state, P: the total pt
// => labels the occupied single particle state of each Fock state
// [Fock_state]: MxR, integer,  M: the total Fock state, R: the total single pt state
// => the Fermion Fock State
// **** VERSION ****
// 1/5/2014 FIRST BUILT 
// 2/24/2014  Modify to accomdate no particle case, ex:[3,0;4,2]
// **** COMMENT ****
// ex:[4,2; 5,3] means 4 states with 2 pts and 5 states with 3 pts.
// index_form='total' (default, label by total), 'site' (label by site)

function [Fock_index,Fock_state]=PIL_Fock_state(occup_num,index_form)
    tot_site=length(occup_num(:,1));
    
    // check input data
    for n=1:tot_site
        if occup_num(n,2) > occup_num(n,1) then
            disp('Error! pt # cannot be larger than state #!');
        end
    end

    [lhs,rhs]=argn(); 
    if rhs==1 then
        index_form='total';
    end
    // Check all hole case
    hole_checker=0;
    counter=0;
    for n=1:tot_site
        if occup_num(n,2)==0 then
            counter=counter+1
            hole_checker(counter)=n;
            occup_num(n,2)=occup_num(n,1);
        end
    end

    // define variables
    tot_comb=1;
    for n=1:tot_site
        tot_comb=tot_comb*PIL_combinatorial(occup_num(n,1),occup_num(n,2));
    end

    // generate dynamic nested for loop commands in text form
    count=0;
    count_m=0;
    loop_begin=[];
    loop_end=[];
    record_index=zeros(tot_comb,sum(occup_num(:,2)));
    for n=1:tot_site
        n0=0;
        for m=1:occup_num(n,2)
            count_m=count_m+1;
            if m==1 then
                loop_begin=loop_begin+' for n'+string(count_m)+'=1:'+string(occup_num(n,1));
            else
                loop_begin=loop_begin+' for n'+string(count_m)+'=n'+string(count_m-1)+'+1:'+string(occup_num(n,1));
            end
            loop_do(count_m)='n'+string(count_m);
            loop_end=loop_end+' end';
        end
        loop_tot=loop_begin+' count=count+1;'+' record_index'+'(count,:)'+'=eval(loop_do'');'+loop_end;
    end
    execstr(loop_tot);

    // recover Fock state form
    Fock_state=zeros(tot_comb,sum(occup_num(:,1)));
    record_index_m=zeros(record_index);
    count_state=0;
    count_pt=0;
    for n=1:tot_site
        record_index_m(:,count_pt+1:count_pt+occup_num(n,2))=record_index(:,count_pt+1:count_pt+occup_num(n,2))+count_state;
        count_pt=count_pt+occup_num(n,2);
        count_state=count_state+occup_num(n,1);
    end

    for n=1:tot_comb
        for m=1:length(record_index_m(1,:))
            Fock_state(n,record_index_m(n,m))=1;
        end
    end

    select index_form
    case 'total' // label by total states
        Fock_index=record_index_m;
    case 'site' // label by each site
        Fock_index=record_index;
    end
    // Change hole output
    if sum(abs(hole_checker))~=0 then
        // Fock_index
        for n=1:length(hole_checker)
            if hole_checker(n)~=1 then
                Fock_index(:,sum(occup_num(1:hole_checker(n)-1,2))+1:sum(occup_num(1:hole_checker(n),2)))...
                =-1*Fock_index(:,sum(occup_num(1:hole_checker(n)-1,2))+1:sum(occup_num(1:hole_checker(n),2)));
                
                Fock_state(:,sum(occup_num(1:hole_checker(n)-1,1))+1:sum(occup_num(1:hole_checker(n),1)))...
                =0*Fock_state(:,sum(occup_num(1:hole_checker(n)-1,1))+1:sum(occup_num(1:hole_checker(n),1)));
            else
                Fock_index(:,1:occup_num(1,2))=-1*Fock_index(:,1:occup_num(1,2));
                Fock_state(:,1:occup_num(1,1))=0*Fock_state(:,1:occup_num(1,1));
            end
        end
        counter=0;
        for n=1:length(Fock_index(1,:))
            if Fock_index(1,n)>0 then
                counter=counter+1;
                Fock_index_tmp(:,counter)=Fock_index(:,n);
            end
        end
        Fock_index=Fock_index_tmp;
    end
endfunction
