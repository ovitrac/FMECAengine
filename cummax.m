function res=cummax(a)
%CUMMMAX cummulative max. cummax is to max what cumsum is to sum.
% ceil(log(size(l(2)))/log(2)) appels recursifs, chaque appel a une complexite de prod(size(a)).
% if coded in java, may be faster with a standard for loop.
% add l(3:end) as third argument of a and of res ; add support for second argument to cummax, using permute of matlab.
% INRA\Daniel Goujot

  l=size(a);
  if l(1)>2 && l(2)>0
    m=ceil(l(1)/2+.7);
    res=cummax([a(1:m-1,:),[a(m:end,:);zeros(mod(l(1),2),l(2))]]);
    res=[res(:,1:l(2));max(res(:,l(2)+1:end),kron(res(end,1:l(2)),ones(m-1,1)))];
    if mod(l(1),2)>0
      res=res(1:end-1,:);
    end
  else
    if l(1)==2
      res=[a(1,:);max(a)];
    else
      if l(2)>1
        res=cummax(a')';
      else
        res=a;
      end
    end
  end
