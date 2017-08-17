// **** Purpose ****
// This code generate the Time-Reversal pairs from ham.state_info 
// **** Variables ****
// [state_info]: 3x3, integer
// <= state_info in ham
// [state_type]: 1x1, string, 'WF', 'SK', 'SKF'
// <= state type of the state_info.
//    'WF': wannier function
//    'SK': Slaster-Koster (must be in cubic harmonic)
//    'SKF': Slaster-Koster Floquet (must be in cubic harmonic)
// [TR_pair]: tot_state/2 x 2, integer
// => each row means a TR pair. 
// **** Version ****
// 04/27/2016 first built
// **** Comment ****

function TR_pair=PIL_TR_pair(state_info, state_type)
    // built SK/SKF suborb flip table 
    function [s_out]=spin_flip(s_in)
        orb_table=[1,2;3,6;4,7;5,8;9,14;10,15;11,16;12,17;13,18;...
        19,26;20,27;21,28;22,29;23,30;24,31;25,32];
        [r,c]=find(orb_table==s_in);
        select c
        case 1
            s_out=orb_table(r,2);
        case 2
            s_out=orb_table(r,1);
        else
            disp('Error: PIL_TR_op, spin_flip does not work!');
        end
    endfunction

    // generate TR_pair
    tot_state=length(state_info(:,1));
    TR_pair=zeros(tot_state/2,2);
    count=0;
    select state_type
    case 'Wannier-Function'
        for n=1:tot_state
            if find(TR_pair(:)==n)==[] then
                count=count+1;
                TR_state=state_info(n,:);
                TR_state($)=-state_info(n,$);
                TR_index=PIL_row_find(state_info,TR_state);
                if TR_state($)==0 then
                    disp('Error: PIL_TR_pair, no spinor!');
                    abort; 
                end
                TR_pair(count,:)=[n,TR_index];
            end
        end
    case 'Slaster-Koster'
        for n=1:tot_state
            if find(TR_pair(:)==n)==[] then
                count=count+1;
                TR_state=state_info(n,:);
                TR_state($)=spin_flip(state_info(n,$));
                TR_index=PIL_row_find(state_info,TR_state);
                if TR_index==[] then
                    disp('Error: PIL_TR_pair, no spinor!');
                    abort; 
                end
                TR_pair(count,:)=[n,TR_index];
            end
        end
    case 'Slaster-Koster Floquet'
        for n=1:tot_state
            if find(TR_pair(:)==n)==[] then
                count=count+1;
                TR_state=state_info(n,:);
                TR_state($)=spin_flip(state_info(n,$));
                TR_index=PIL_row_find(state_info,TR_state);
                if TR_index==[] then
                    disp('Error: PIL_TR_pair, no spinor!');
                    abort; 
                end
                TR_pair(count,:)=[n,TR_index];
            end
        end
    else
        disp('Error: PIL_TR_pair, state_type is wrong!');
        abort
    end
endfunction
