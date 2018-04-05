function [ i,v] = max_such_that(A,f)
    N=length(A);
    for i=1:N
        if ~f(A(i))
            i=i-1;
            break
        end
    end
    if i<1
        v=[];
    else
        v=A(i);
    end
end

