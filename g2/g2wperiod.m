function w = g2wperiod(w)
% function w = g2wperiod(w)
% fix periodicity for extra cells at left/right in phi (2nd index)
w(:,  1,:) = w(:,end-1,:);
w(:,end,:) = w(:,    2,:);
%EOF
