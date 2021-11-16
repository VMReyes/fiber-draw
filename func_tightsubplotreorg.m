function []=  func_tightsubplotreorg(myfighandle,Nh, Nw)

%if(~isempty(myfighandle))
%    figure(myfighandle)
%end

gap = .02*2;
 marg_h = .05; 
 marg_w = .05; 
 
 gap = .02*2;
 marg_h = .08; 
 marg_w = .08;
 
gap = [gap gap];
marg_w = [marg_w marg_w];
 marg_h = [marg_h marg_h];

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
py = 1-marg_h(2)-axh; 



ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        mysubfigHANDLE0 = subplot(Nh,Nw,ii);
        set(mysubfigHANDLE0,'Position',[px py axw axh])
%        set(allAxesInFigure(ii),'Position',[px py axw axh])
%         ha(ii) = axes('Units','normalized', ...
%             'Position',[px py axw axh], ...
%             'XTickLabel','', ...
%             'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end




