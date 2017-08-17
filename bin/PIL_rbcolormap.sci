// this function generates colormap from blue to red
function rbcolormap=PIL_rbcolormap(tot_color)
    rbcolormap=PIL_k_path([1 0 0; 0 0 1],tot_color);
endfunction
