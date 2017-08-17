// **** Purpose ****
// convert low-dimension vectors to three dimension representation
// **** Variables ****
// [vec_list_in]: nx3 | nx2 | nx1, real
// <= r_points in row
// [vec_list_out]: nx3, real
// => three dimension r_points
// **** Version ****
// 05/29/2014
// **** Comment ****
function vec_list_out=PIL_vec_3d(vec_list_in)
    tot_vec=length(vec_list_in(:,1));
    vec_dim=length(vec_list_in(1,:));
    select vec_dim
    case 1
        vec_list_out=cat(2,vec_list_in,zeros(tot_vec,2));
    case 2
        vec_list_out=cat(2,vec_list_in,zeros(tot_vec,1));
    case 3
        vec_list_out=vec_list_in;
    end
endfunction
