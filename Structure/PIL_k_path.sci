// **** Purpose ****
// generates the k-point along your assigned k-path.
// **** Variables ****
// [k_path]: nx3, real
// <= k points to specify your path. Can be many points.
// [k_div]: 1x1, integer
// <= how many divides for each path
// [lat_const]: 1x1, real
// <= the lattice constant, default=1
// [k_div_opt]: 'unit' or 'seg'
// <= how to divide each k-path, default='seg'
// [k_path_point]: N x 1, integer
// => the number of divisions of each k-path
// [k_path_point]: k_div x 3, real
// => the k_path ponts
// **** Version ****
// 05/01/2014 1st version
// 10/20/2014 add k_div_opt function
// 04/16/2014 delete duplicated k-points
// **** Comment ****
// this function provides you two difference way to generate k-points along
// a k-path. if you chose 'unit', you should input lat_const to determine
// the unit length in k-space. So any BZ has 2%pi as maxima and independent
// of lattice constant.  
function [k_path_point,k_path_div]=PIL_k_path(k_path,k_div,k_div_opt,lat_const) 
    [lhs,rhs]=argn();
    select rhs
    case 2 
        k_div_opt='seg';
        lat_const=1;
    case 3
        lat_const=1;
    end

    tot_seg=length(k_path(:,1))-1;
    select k_div_opt
    case 'seg'
        k_path_div=k_div*ones(tot_seg,1);
    case 'unit'
        for n=1:tot_seg
            k_path_div(n)=..
            round(k_div*norm(lat_const*(k_path(n+1,:)-k_path(n,:))));
        end
    end

    k_path_point=zeros(sum(k_path_div)-length(k_path_div)+1,3);
    count=0;
    for n=1:tot_seg
        // prevent k-point double counting
        if n==1 then
            for m=1:k_path_div(n)
                count=count+1;
                k_path_point(count,:)=k_path(n,:)..
                +((m-1)/(k_path_div(n)-1))*(k_path(n+1,:)-k_path(n,:));
            end
        else
            for m=2:k_path_div(n)
                count=count+1;
                k_path_point(count,:)=k_path(n,:)..
                +((m-1)/(k_path_div(n)-1))*(k_path(n+1,:)-k_path(n,:));
            end
        end
    end
    // fix k_path_div double counting
    k_path_div(2:$)=k_path_div(2:$)-1;
endfunction 

//example:
//k_path=[2*%pi/a 0 0;0 0 0;0 0 2*%pi/a];
//k_div=100;
//PIL_k_path(k_path,k_div) 
