// **** Purpose ****
// generate the index of a dynamic nested for loops
// **** Variables ****
// [loop_mat]: n x 2, integer
// <= describes the range of each loop
//    ex:[2,5;3,7;1,5]  means 2:5, 3:7, 1:5
// [loop_index]: n x m, integer
// => the index of each loop
// **** Version ****
// 05/01/2014 first built
// 05/20/2014 fix bug, when length(loop_mat)==1 has error!
// 05/30/2014 fullly rewrite to accept all kinds of dynamical nested loops
// 01/24/2016 fully rewrite using matrix reshape algorithm. super fast now!
// **** Comment ****
function loop_index=PIL_nest_loop(loop_mat)
    tot_loop=length(loop_mat(:,1));
    loop_new=flipdim(loop_mat(:,2)-loop_mat(:,1)+1,1);
    loop_index=zeros(prod(loop_new),tot_loop);
    for n=1:tot_loop
        if n==1 then
            tmp_mat=repmat(1:loop_new(n),1,prod(loop_new(n+1:$)));
        elseif n==tot_loop
            tmp_mat=repmat(1:loop_new(n),prod(loop_new(1:n-1)),1);
        else
            tmp_mat=repmat(1:loop_new(n),prod(loop_new(1:n-1)),prod(loop_new(n+1:$)));
        end
        loop_index(:,n)=tmp_mat(:);
        loop_index(:,n)=loop_index(:,n)-1+loop_mat(tot_loop-n+1,1)
    end
    loop_index=flipdim(loop_index,2);
endfunction
