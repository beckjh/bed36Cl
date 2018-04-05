function results2 = thin(results,factor)
results2 = results;
results2.state_container = containers.Map('keyType','uint64','valueType','any');
for n = 1:length(results.repeat_history)
    results2.repeat_history{n} = results.repeat_history{n}(1:factor:end);
    results2.state_history{n} = results.state_history{n}(1:factor:end);
    for i = 1:length(results2.state_history{n})
        id = results2.state_history{n}{i};
        results2.state_container(id) = results.state_container(id);
    end
end
end
