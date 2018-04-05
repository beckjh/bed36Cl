function y = tensor_evaluation(values,tensor,tensor_bounds)
    indices = value2index(values,tensor_bounds);
    y = tensor(sub2ind(size(tensor),indices(:,1),indices(:,2),indices(:,3),indices(:,4)));
    y = y(:);
end
function i = value2index(values,value_bounds)
    i = zeros(size(values));
    for d = 1:size(value_bounds,1)
        len = value_bounds(d,3);
        i(:,d) = max(1,min(round((values(:,d)-value_bounds(d,1))/(value_bounds(d,2)-value_bounds(d,1))*len),len));   
    end
end