figHandles = findobj('Type', 'figure');

for ff = 1:length(figHandles)
    set(figHandles(ff),'units','normalized','outerposition',[0 0 1 1])

    %tit = ['ttfigs' num2str(figHandles(ff))];
    tit = ['ttfigs' num2str(figHandles(ff).Number)];
    %print(figHandles(ff), tit,'-dtiffnocompression' )
    print(figHandles(ff), tit,'-dtiffn' )
end