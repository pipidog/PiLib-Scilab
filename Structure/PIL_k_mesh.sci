// **** Purpose ****
// Give four points in the k-space. The first one is origin point and
// the other three points span a space in the k-space. This code will 
// perform k-mesh on this hexahedral. 
// **** Variables ****
// [rec_vec]: 3x3, real, if [], then eye(3,3)
// <= reciprocal row vectros. If [], then eye(3,3). It is used if span_vec
//    are in reduced coordinate the you want to k-points to output in
//    reduced coordinate.
// [k_mesh]: 1x3, int
// <= mesh along each axis. Note: if any k_mesh=1, there will be no
//    mesh on that axis. If any k_mesh=2, there will be only origin and
//    end point in that direction.
// [span_vec]: 4x3, real, default: [0,0,0;eye(3,3)]
// <= give 4 points to define a hexahedral in k-space. The first one
//    will be the origin. The other three span the space. 
//    If [], will mesh the whole BZ.
// [k_format]: 1x1, string, 'red' / 'cart', default: 'red'
// <= define the meaning of the span_vec variable. reduced or cartisian
//    coordinate. 
// [end_point]: 1x1, string , 'on' / 'off', default:'off'
// <= whether to include the end_point. Sometimes the end_point are 
//    equivalent to the origin, e.g 0*pi/a and 2*pi/a. Therefore, one
//    can choose whether to include the end_point. If 'off', the mesh 
//    will not include end_point.  
// [k_point]: nx3, real
// => the k-points. 
// **** Version ****
// 05/01/2014 first version
// 06/03/2014 fully rewrite, more concise, accept low dimensions
// 05/02/2016 fully rewrite to perform more general k-mesh. 
// **** Comment ****
// If one just want to perform a general k-mesh on the whole BZ, then
// just use PIL_k_mesh(rec_vec,k_mesh).

function [k_point]=PIL_k_mesh(rec_vec,k_mesh,span_vec,k_format,end_point)
    [lhr,rhs]=argn()
    select rhs
    case 2
        span_vec=[0,0,0;eye(3,3)]
        k_format='red';
        end_point='off';
    case 3
        k_format='red';
        end_point='off';
    case 4
        end_point='off';
    end
    
    // construct difference vectors 
    diff_vec=zeros(3,3);
    k_div=PIL_nest_loop(cat(2,[0,0,0]',(k_mesh-1)'));
    for n=1:3
        diff_vec(n,:)=span_vec(n+1,:)-span_vec(1,:);
        if k_mesh(n)~=1 then
            select end_point
            case 'on'
                k_div(:,n)=k_div(:,n)/(k_mesh(n)-1);    
            case 'off'
                k_div(:,n)=k_div(:,n)/(k_mesh(n));
            end    
        end
    end

    // construct k-points
    k_point=repmat(span_vec(1,:),length(k_div(:,1)),1)+k_div*diff_vec;
    if k_format=='red' then
        k_point=k_point*rec_vec;
    end
endfunction
