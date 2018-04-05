function state=strip_pre_erosion(state)
    i_eroded=find(state.times>=state.T_init,1);
    state.times = state.times(i_eroded:end);
    state.jumps = state.jumps(i_eroded:end);
end