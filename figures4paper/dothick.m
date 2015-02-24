% dothick.m
% set thick lines
% requires: thick = 1; 
%  xlab for xlabel    'text'
%  ylab for ylabel
%  titl for title
%  legn for legend {'first','second'}, legloc for location
%  colb for colorbar 'vert'
% _____________________________________________________________________________
do_thk = 0; if(exist('thick')), if(thick ~= 0), do_thk = 1; end; end
if(do_thk)
  if(thick > 1)       lWid = 3; lMar = 16;
  elseif(thick > 0.8) lWid = 2; lMar = 12;
  else                lWid = 1; lMar = 7;end 
  set(findobj('type', 'line'), 'LineWidth',lWid,'MarkerSize',lMar);
  if(exist('nocol') & nocol == 1)
    fprintf(1,'nocol ..\n');
    set(gca,...
      'FontSize',16,'FontWeight','normal','FontName','times',...
      'LineWidth',lWid);
  else
    set(gca,'ColorOrder',[0 0.4 0;0.4 0 0;0 0 0.4;0.4 0.4 0;0 0.4 0.4;0.4 0 0.4;0 0 0],...
      'FontSize',16,'FontWeight','normal','FontName','times',...
      'LineWidth',lWid);
  end
  if(exist('xlab') & ~isempty(xlab))
    xlabel(xlab,'FontSize',16,'FontWeight','normal','FontName','times');
  end
  if(exist('ylab') & ~isempty(ylab)) 
    ylabel(ylab,'FontSize',16,'FontWeight','normal','FontName','times');
  end
  if(exist('titl') & ~isempty(titl)) 
     title(titl,'FontSize',16,'FontWeight','normal','FontName','times'); 
  end
  if(exist('colb') & ~isempty(colb)) 
     colorbar(colb,'LineWidth',lWid,'FontSize',16,...
       'FontWeight','normal','FontName','times'); 
  end
  if(exist('legn') & ~isempty(legn)),
    mylegloc = 'Best';
    if(exist('legloc') & ~isempty(legloc)), mylegloc = legloc; end
    if(exist('leghandle') & ~isempty(leghandle)),
      % leghandle
      legend(leghandle,legn,'FontSize',16,'FontWeight','normal',...
		'FontName','times','Location',mylegloc); %
    else
      leghand2 = legend(legn,'FontSize',16,'FontWeight','normal',...
		'FontName','times','Location',mylegloc); % ,'Color','none'); %
      %,'interpreter','latex');
      if(0 & exist('legcol','var') & ~isempty(legcol)),
        legend(leghand2,'Color',legcol);
      end
    end
  end
else % ________________________________________________________________________
  if(exist('xlab') & ~isempty(xlab)) xlabel(xlab); end   
  if(exist('ylab') & ~isempty(ylab)) ylabel(ylab); end   
  if(exist('titl') & ~isempty(titl)) title(titl);  end
  if(exist('colb') & ~isempty(colb)) colorbar(colb); end
  if(exist('legn') & ~isempty(legn)) 
    mylegloc = 'Best';
    if(exist('legloc') & ~isempty(legloc)), mylegloc = legloc; end
    if(exist('leghandle') & ~isempty(leghandle)),
       legend(leghandle,legn,'Location',mylegloc);
    else legend(legn,'Location',mylegloc); end
  end
end  
% reset
if( exist('no_reset') & no_reset > 0 )
else xlab = []; ylab = []; titl = []; legn = []; leghandle = []; colb = []; 
     legloc = [];
end
%
%%EOF
