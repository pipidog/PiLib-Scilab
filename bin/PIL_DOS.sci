// **** Purpose ****
// It counts the energy levels for a specific energy interval, so you 
// can plot DOS. 
// **** variables ****
// [E_level]: nx1, real
// <= the energy levels
// [E_interval]:1x1, real
// <= the energy interval
// [E_point]: nx1, real
// => the energy points
// [E_DOS]: nx1, integer
// => The DOS of each interval
// **** Version ****
// 06/04/2014 1st version
// **** Comment ****
// 1. To draw DOS, just plot(E_point,E_DOS)
// 2. Because only the relative strength are important 
//    for DOS, so it renormalizes the maximum to 1 
function [E_point,E_DOS]=PIL_DOS(E_level,E_interval)
    E_level=gsort(E_level,'g','i');
    tot_div=ceil((E_level($)-E_level(1))/E_interval)+2;
    E_DOS=zeros(1,tot_div);
    E_point=zeros(1,tot_div);
    for n=1:tot_div+2
        E_lower=E_level(1)+(n-2)*E_interval;
        E_upper=E_level(1)+(n-1)*E_interval;
        E_point(n)=(E_lower+E_upper)/2;
        E_DOS(n)=length(find(E_level<E_upper & E_level>E_lower));
    end
    // renormalization
    E_DOS=E_DOS/max(E_DOS);
    title('DOS','fontsize',5);
endfunction
