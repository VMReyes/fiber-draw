function clrstr = func_clrstr(str,n)

if(0)

    str{1} = 'r';
    str{2} = 'b';
    str{3} = 'g';

    %numstr = length(str);
    
    %n = 1:5;
    %nWhichStr = mod(n+numstr-1,numstr)+1;
    n = 1:5;
    clrstr = func_clrstr(str,n)
end

    numstr = length(str);
    clrstr = mod(n+numstr-1,numstr)+1;
    