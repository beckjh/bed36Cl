function [states,repeats] = burn_in(results,p_burn_in,range)
if nargin < 3
    range = 1:results.settings.group_size;
end
state_history = results.state_history(range);
repeat_history = results.repeat_history(range);
state_container = results.state_container;
N_chains = length(state_history);
states = cell(0,1);
repeats = zeros(0,1);
for n = 1:N_chains
   N_states = sum(repeat_history{n});
   burned = 0;
   N_burn = floor(p_burn_in*N_states);
   i = 1;
   i_max = length(repeat_history{n})+1-find(repeat_history{n}(end:-1:1),1);
   while burned < N_burn && i < length(repeat_history{n})
      n_burn = min(N_burn-burned, repeat_history{n}(i));
      burned = burned+n_burn;
      if n_burn < repeat_history{n}(i)
        repeat_history{n}(i) = repeat_history{n}(i)-n_burn;
      end
      i = i+1;
   end
   state_history{n} = state_history{n}(i:i_max);
   repeat_history{n} = repeat_history{n}(i:i_max);
   state_history_mat = cell2mat(state_history{n});
   [state_ids,~,u] = unique(state_history_mat);
   new_repeats = zeros(length(state_ids),1);
   new_states = cell(length(state_ids),1);
   for i_id = 1:length(state_ids)
       new_states{i_id} = state_container(state_ids(i_id));
   end
   for i = 1:length(state_history_mat)
       new_repeats(u(i)) = new_repeats(u(i))+repeat_history{n}(i);
   end
   states = cat(1,states,new_states);
   repeats = cat(1,repeats,new_repeats);
end
end