function ndif = g2pardif(gparold,gparnew)
% function ndif = g2pardif(gparold,gparnew)
% prints differences between old and new parameters in structs gparold,gparnew
names0 = fieldnames(gparold); ndif = 0;
for k=1:length(names0),
  namek = char(names0(k)); gv0 = gparold.(namek); gv1 = gparnew.(namek); 
  if(isnumeric(gv0)),
    if(gv0 ~= gv1), 
      fprintf(1,'  %s: %g -> %g\n',namek,gv0,gv1); ndif = ndif+1;
    end
  elseif(ischar(gv0)),
    if(~strcmp(gv0,gv1)),
      fprintf(1,'  %s: %s -> %s\n',namek,gv0,gv1); ndif = ndif+1;
    end
  else fprintf(1,'g2pardif.m -- unexpected type\n');
  end
end
%EOF
