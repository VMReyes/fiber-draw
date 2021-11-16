function [dataaout,timeout] = func_simplefilter(dataa,fltLEN, SHAPE)

if(nargin < 3)
    SHAPE = 'valid';
end

%fltLEN - odd

[m,n] = size(dataa);

hlflen = (fltLEN - 1)/2;
B = ones(1, fltLEN);
%B =  window(@gausswin,fltLEN);
B = B/sum(sqrt(B.^2));

%filtering (and then, seperately, scaling)
if(1)
    if(strcmp(SHAPE,'valid'))
        dataaTEMP = zeros(length((hlflen+1):(m-hlflen)),n);
    end
    if(strcmp(SHAPE,'same'))
        dataaTEMP = zeros(m,n);
    end
    %dataaSCALED = dataa((hlflen+1):(end-hlflen),:);
    %clear theSCALEsave;
    for nDataCol = 1:n
        %perform filter
%        flted = conv(dataa(:,nDataCol),B, 'valid');
        flted = conv(dataa(:,nDataCol),B, SHAPE);
        dataaTEMP(:,nDataCol) = flted;
    end
    dataaout = dataaTEMP;
    sz = size(dataaout);
    len = sz(1);
     
    if(strcmp(SHAPE,'valid'))
        timeout = (hlflen+1):(m-hlflen);
    end     
    if(strcmp(SHAPE,'same'))
        timeout = 1:m;
    end
    
end