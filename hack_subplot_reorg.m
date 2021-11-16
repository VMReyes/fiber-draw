% load testsavefig50000.mat
% allAxesInFigure = findall(myfigHANDLE,'type','axes');
% 
% %see https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
% 
% Nh = 4;
% Nw = 5;
% 
% gap = .02*2;
%  marg_h = .05; 
%  marg_w = .05; 
% gap = [gap gap];
% marg_w = [marg_w marg_w];
%  marg_h = [marg_h marg_h];
% 
% axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
% axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
% py = 1-marg_h(2)-axh; 
% 
% 
% 
% ii = 0;
% for ih = 1:Nh
%     px = marg_w(1);
%     
%     for ix = 1:Nw
%         ii = ii+1;
%         mysubfigHANDLE0 = subplot(Nh,Nw,ii);
%         set(mysubfigHANDLE0,'Position',[px py axw axh])
% %        set(allAxesInFigure(ii),'Position',[px py axw axh])
% %         ha(ii) = axes('Units','normalized', ...
% %             'Position',[px py axw axh], ...
% %             'XTickLabel','', ...
% %             'YTickLabel','');
%         px = px+axw+gap(2);
%     end
%     py = py-axh-gap(1);
% end
% 
% 
% %%
% 
% load testsavefigs.mat;

%%
% for numfigHs = 1:length(myfigHANDLEarray)
%     
%     myfigHANDLE = myfigHANDLEarray{numfigHs};
%     allAxesInFigure = findall(myfigHANDLE,'type','axes');
%     nsubs = length(allAxesInFigure)/4;
%     
%     func_tightsubplotreorg(myfigHANDLE,4, nsubs);
% end

%%


%%

for numfigHs = 1:length(myfigHANDLEarray)
 
    %newfig = figure(myfigNUMBERarray(numfigHs));
    myfigHANDLE = myfigHANDLEarray{numfigHs};
    figure(myfigHANDLE)
    
%     allAxesInFigure = findall(myfigHANDLE,'type','axes');
%     nsubs = length(allAxesInFigure)/4;
% 
%     func_tightsubplotreorg(myfigHANDLE,4, nsubs);
    
    submn = myfigSUBPLOTarray(:,numfigHs);
    func_tightsubplotreorg(myfigHANDLE,submn(1), submn(2));

end


%%
for numfigHs = 1:length(myfigHANDLEarray)
    myfigHANDLE = myfigHANDLEarray{numfigHs};
    newfig = figure(myfigNUMBERarray(numfigHs));
    copyobj(allchild(myfigHANDLE),newfig)
    
     
    %allAxesInFigure = findall(myfigHANDLE,'type','axes');
    %nsubs = length(allAxesInFigure)/4;

    %func_tightsubplotreorg(newfig,4, nsubs);
end

%%
for numfigHs = 1:length(myfigHANDLEarray)
 
    %newfig = figure(myfigNUMBERarray(numfigHs));
    myfigHANDLE = myfigHANDLEarray{numfigHs};
    close(myfigHANDLE)

end

%%


% ii = 0;
% for ih = 1:nsubs
%     px = marg_w(1);
%     
%     for ix = 1:Nw
%         ii = ii+1;
%         set(allAxesInFigure(ii),'Position',[px py axw axh])
% %         ha(ii) = axes('Units','normalized', ...
% %             'Position',[px py axw axh], ...
% %             'XTickLabel','', ...
% %             'YTickLabel','');
%         px = px+axw+gap(2);
%     end
%     py = py-axh-gap(1);
% end
