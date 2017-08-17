// **** Purpose ****
// print a data array in a fixed format
// **** Variables ****
// desc: string
// <= a string describes this data
// [A]: nxn, real or complex
// <= the data array to be print out
// [A_type]: char, 'r' (real), 'i' (integer), 'c' (complex), 'sp' (sparse)
// <= specify the format of data array A
// [fid]: 1x1, integer
// <= the file ID of the print out file 
// [prt_order]: 1x1, 'on' / 'off', optional, default='on'
// <= if 'on', use 0.? * 10^n representation 
//    if 'off', use normal representation. 
//    'off' is particular useful for variables that combine integers and
//    float numbers. 
// **** Version ****
// 05/01/2014 1st version
// 06/05/2014 add line-text print, add data_type output
// 06/07/2014 add sparse matrix format
// 06/07/2014 add text matrix format (separated by '#')
// 04/26/2015 print row numbers for 'i', 'r', 'c' cases
// 02/17/2016 add prt_order to control order representation
// 04/16/2016 add 'sp' to accept native scilab sparse constant matrix
// **** Comment ****
// This fnction exports your variable using the standard formate of my Fortran
// PI library. desc is a string description that appears in your output file.
// A is the output matrix. A_type can be 'r'(real float), 'i'(integer),'c'(complex)
// that specify your output format. fid is your file ID. 

function PIL_print_mat(desc,A,A_type,fid,prt_order,py_format)
    [lhs,rhs]=argn()
    select rhs
    case 4
        prt_order='on'
        py_format='off'
    case 5
        py_format='off'
    end
    if py_format~='on' & py_format~='off' then
        py_format='off'
    end

    mputl('',fid);
    select A_type
    case 'i'
        A_dim=size(A)
        order=0
        // add row number to A
        A=A*10^(-order);
        A=cat(2,[1:A_dim(1)]',A);

        out_format=strcat(repmat('  %6d',1,A_dim(2)))
        column_format=strcat(repmat('  %6d',1,A_dim(2)))
        mfprintf(fid,'============= PiLib Variable =============\n')
        mfprintf(fid,'%s\n',desc)
        mfprintf(fid,'ORDER= %4d, SIZE=[%6d,%6d], TYPE=%s\n',..
        order,A_dim(1),A_dim(2),'INTEGER')
        mfprintf(fid,'\n')
        mfprintf(fid,'      '+column_format+'\n',1:A_dim(2))
        if py_format=='off' then
            mfprintf(fid,'\n')
        end
        mfprintf(fid,'%6d'+out_format+'\n',A);
    case 'r'
        A_dim=size(A)
        select prt_order
        case 'on'
            order=max(abs(A));
            if order~=0 then
                order=floor(log10(order));
            end
        case 'off'
            order=0 
        end

        // add row number to A
        A=A*10^(-order);
        A=cat(2,[1:A_dim(1)]',A);

        out_format=strcat(repmat('  %10.6f',1,A_dim(2)))
        column_format=strcat(repmat('  %10d',1,A_dim(2)))
        mfprintf(fid,'============= PiLib Variable =============\n')
        mfprintf(fid,'%s\n',desc)
        mfprintf(fid,'ORDER= %4d, SIZE=[%6d,%6d], TYPE=%s\n',..
        order,A_dim(1),A_dim(2),'REAL')
        mfprintf(fid,'\n')
        mfprintf(fid,'      '+column_format+'\n',1:A_dim(2))
        if py_format=='off' then
            mfprintf(fid,'\n')
        end
        mfprintf(fid,'%6d'+out_format+'\n',A);
    case 'c'
        A_dim=size(A)
        select prt_order
        case 'on'
            order=max(abs(A));
            if order~=0 then
                order=floor(log10(order));
            end
        case 'off'
            order=0 
        end
        // sort to A to two column real-complex form
        B=zeros(A_dim(1),2*A_dim(2));
        B(:,[1:2:$])=real(A);
        B(:,[2:2:$])=imag(A);
        A=B;

        // add row number to A
        A=A*10^(-order);
        A=cat(2,[1:A_dim(1)]',A);


        column_format=strcat(repmat('%20d',1,A_dim(2)));
        out_format=strcat(repmat('%10.6f',1,2*A_dim(2)));

        mfprintf(fid,'============= PiLib Variable =============\n')
        mfprintf(fid,'%s\n',desc)
        mfprintf(fid,'ORDER= %4d, SIZE=[%6d,%6d], TYPE=%s\n',..
        order,A_dim(1),A_dim(2),'COMPLEX')
        mfprintf(fid,'\n')
        mfprintf(fid,'        '+column_format+'\n',1:A_dim(2))
        if py_format=='off' then
            mfprintf(fid,'\n')
        end
        mfprintf(fid,'%6d  '+out_format+'\n',A);
    case 'sp'
        [ij,v,mn]=spget(A);
        A=cat(1,[mn,0],cat(2,ij,v))

        A_dim=size(A)
        select prt_order
        case 'on'
            order=max(abs(A(:,3)));
            if order~=0 then
                order=floor(log10(order));
            end 
        case 'off'
            order=0 
        end
        column_format='  %6d  %6d  %20d'

        mfprintf(fid,'============= PiLib Variable =============\n')
        mfprintf(fid,'%s\n',desc)
        mfprintf(fid,'ORDER= %4d, SIZE=[%6d,%6d], TYPE=%s\n',..
        order,A_dim(1),A_dim(2),'SPARSE')
        mfprintf(fid,'\n')
        mfprintf(fid,column_format+'\n',1:A_dim(2))
        if py_format=='off' then
            mfprintf(fid,'\n')
        end

        A=cat(2,A,A(:,3));
        A(:,1:2)=round(A(:,1:2));
        A(:,3)=real(A(:,3))*10^(-order);
        A(:,4)=imag(A(:,4))*10^(-order);
        mfprintf(fid,'  %6d  %6d  %10.6f%10.6f\n',A);

    case 's'
        A_dim=size(A);
        order=0
        out_format=strcat(repmat('%s # ',1,A_dim(2)))
        column_format=''
        mfprintf(fid,'============= PiLib Variable =============\n')
        mfprintf(fid,'%s\n',desc)
        mfprintf(fid,'ORDER= %4d, SIZE=[%6d,%6d], TYPE=%s\n',..
        order,A_dim(1),A_dim(2),'STRING')
        mfprintf(fid,'\n')
        if py_format=='on' then
            mfprintf(fid,'\n')
        end
        mfprintf(fid,out_format+'\n',A);
    end
endfunction
