function [ i,v] = min_such_that(A,f)
    N=length(A);
    for i=N:-1:1
        if ~f(A(i))
            i=i+1;
            break
        end
    end
    if i>N
        v=[];
    else
        v=A(i);
    end
end
