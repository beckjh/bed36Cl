function multi_idx = multiidx_gen(L,rule,w,base,multiidx,multi_idx)

if nargin==3
      base = 0;
      multiidx=[];
      multi_idx=[];
elseif nargin==4
      multiidx=[];
      multi_idx=[];
end

if length(multiidx)~=L   
      i=base;
      if isscalar(w)
          while rule([multiidx, i]) <= w
                multi_idx = multiidx_gen(L,rule,w,base,[multiidx, i],multi_idx);
                i=i+1;
          end
      else
          tocontinue=1;
          while tocontinue
                for d=1:1:(length(multiidx)+1)
                    if rule([multiidx, i],d) > w(d)
                        tocontinue=0;
                    end
                end
                if tocontinue
                    multi_idx = multiidx_gen(L,rule,w,base,[multiidx, i],multi_idx);
                    i=i+1;
                end
          end
      end
else  
      multi_idx=[multi_idx; multiidx];
end