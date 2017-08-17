// **** Purpose ****
// This plot the bar diagram of an input coulmn data
// **** Variables ****
// [A]: real, Nx1
// <= data to plot  
// [tot_div]: int, 1x1
// <= total division of the plot
// [data_plot]: 'on' / 'off'
// <= whether to plot the bar diagram
// [bar_data]: real, tot_div x 1
// => data to plot bar diagram
// **** Version ****
// 08/08/2015
function bar_data=PIL_bar_plot(A, tot_div, data_plot)
    data_mesh=linspace(min(A),max(A),tot_div+1)
    bar_data=zeros(tot_div,2);
    for n=1:tot_div
        bar_data(n,1)=data_mesh(n);
        bar_data(n,2)=length(find(A >=data_mesh(n)..
        & A < data_mesh(n+1)))
    end
    if data_plot=='on'
        bar(bar_data(:,1),bar_data(:,2))
    end
endfunction

