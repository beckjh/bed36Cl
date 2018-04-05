classdef Sampler
    properties
       T_min;
       T_max;
       T_total;
       J_total;
       N_bm;
       lambda_interval_bm=4;
       Lambda_2=4;
       state;
       n_proposals;
       prior;
       S;
       no_more_slips;
       d_min;
       d_max;
       parameters
       var_parameters;
       cutoff;
    end
    methods   
        function obj=Sampler(fault,state,settings)
           obj.cutoff=settings.cutoff;
           obj.prior=fault.prior;
           obj.T_min=settings.T_min;
           obj.T_max=settings.T_max;
           obj.T_total=obj.T_max-obj.T_min;
           obj.parameters=fault.parameters;
           obj.J_total=fault.parameters.H_sc;
           obj.N_bm=2^(ceil(log(obj.T_total/settings.dT)/log(2)))+1;
           obj.no_more_slips=fault.prior.no_more_slips;
           obj.d_min=fault.prior.d_min;
           obj.d_max=fault.prior.d_max;
           if isempty(state)
               obj.state.S=unifrnd(0,1);
               obj.state.tau=(obj.prior.tau_max+obj.prior.tau_min)/2;
               MH_factor=0;
               obj.state.p_zero_freq=unifrnd(0,obj.prior.p_zero_freq_max);
               while MH_factor==0
                   obj.state.lambda_switches=unifrnd(obj.prior.lambda_switches_min,obj.prior.lambda_switches_max);
                   N_switches=poissrnd(obj.state.lambda_switches);
                   obj.state.freq_switches=sort(unifrnd(obj.T_min,obj.T_max,[N_switches,1]));
                   obj.state.freq_values=zeros(N_switches+1,1);
                   obj.state.freq_mean = unifrnd(obj.prior.freq_mean_min,obj.prior.freq_mean_max);
                   obj.state.freq_alpha = unifrnd(obj.prior.freq_alpha_min,obj.prior.freq_alpha_max);
                   for j=1:N_switches+1
                       if unifrnd(0,1)>=obj.state.p_zero_freq
                           obj.state.freq_values(j) = invgamrnd(obj.state.freq_mean,obj.state.freq_alpha);
                       end
                   end
                   obj.state.T_init=unifrnd(obj.prior.T_init_min,obj.prior.T_init_max);
                   obj.state.bm=zeros(obj.N_bm,1);
                   obj.state.bm=obj.generate_bm(obj.state.bm,1);
                   [obj.state,MH_factor]=obj.add_slip_history(obj.state,obj.state,1);
               end           
               assert(MH_factor>0)
           else
               obj.state=state;
           end
           obj.var_parameters=[];
           parameter_names=fieldnames(obj.parameters);
           for i=1:length(parameter_names)
                if isa(obj.parameters.(parameter_names{i}),'function_handle')
                    obj.var_parameters=[obj.var_parameters;i];
                    if isempty(state)
                        obj.state.(parameter_names{i})=obj.parameters.(parameter_names{i})([]);
                    end
                else
                   obj.state.(parameter_names{i})=obj.parameters.(parameter_names{i});
                end
           end
           assert(obj.is_feasible_freq(obj.state.freq_values,obj.state));
           assert(obj.is_feasible_slip_history(obj.state));
           assert(length(obj.state.times)==length(obj.state.jumps));
           obj.n_proposals=28;
        end
        function [proposal,MH_ratio,type]=propose(obj)
            proposal=obj.state;
            MH_ratio=1;
            n_combine=geornd(0.5)+1;
            for combi=1:n_combine
                if combi==1
                    weights=[1,1,2,1,4,1,1,1,4,1,1,1,1,1];
                else
                    weights=ones(1,obj.n_proposals/2);
                end
                group=discrete_draw(weights,1);
                subtype=unidrnd(2);%global/local
                type=2*(group-1)+subtype;
                switch group
                    case 1%1 2
                        [proposal,MH_factor]=propose_jump_sizes(obj,subtype,proposal);
                    case 2%3 4
                        [proposal,MH_factor]=propose_bm(obj,subtype,proposal);
                    case 3%5 6
                        [proposal,MH_factor]=propose_freq_switches(obj,subtype,proposal);
                    case 4%7 8
                        [proposal,MH_factor]=propose_freq_values(obj,subtype,proposal);
                    case 5%9 10    
                        [proposal,MH_factor]=propose_number_freq_switches(obj,subtype,proposal);    
                    case 6%11 12
                        [proposal,MH_factor]=propose_tau(obj,subtype,proposal);
                    case 7%13 14
                        [proposal,MH_factor]=propose_S(obj,subtype,proposal);
                    case 8% 15 16                        
                         [proposal,MH_factor]=propose_T_init(obj,subtype,proposal);
                    case 9% 17 18
                         [proposal,MH_factor]=propose_freq_interlude(obj,subtype,proposal);  
                    case 10% 19 20
                         [proposal,MH_factor]=propose_lambda_switches(obj,subtype,proposal);
                    case 11% 21 22
                        [proposal,MH_factor]=propose_parameters(obj,subtype,proposal);
                    case 12%23 24
                        [proposal,MH_factor]=propose_freq_mean(obj,subtype,proposal);
                    case 13%25 26
                        [proposal,MH_factor] = propose_freq_alpha(obj,subtype,proposal);
                    case 14%27 28
                        [proposal,MH_factor] = propose_p_zero_freq(obj,subtype,proposal);
                end
                if max(proposal.times)>obj.no_more_slips
                    MH_factor=0;
                end
                MH_ratio=MH_ratio*MH_factor;
                if isnan(MH_ratio)
                    MH_ratio=0;
                end
                if MH_ratio==0
                    break
                end               
                assert((MH_ratio==0)||obj.is_feasible_freq(proposal.freq_values,proposal));
                assert((MH_ratio==0)||~any(isnan(proposal.jumps)))
                assert((MH_ratio==0)||(length(proposal.times)==length(proposal.jumps)));
            end
            assert((MH_ratio==0)||obj.is_feasible_freq(proposal.freq_values,proposal));
            assert((MH_ratio==0)||~any(isnan(proposal.jumps)))
            assert((MH_ratio==0)||(length(proposal.times)==length(proposal.jumps)));
        end
    end 
    methods
        function [proposal,MH_ratio]=propose_p_zero_freq(obj,subtype,proposal)
           old=proposal;
           if obj.prior.p_zero_freq_max>0
               switch subtype
                   case 1
                       proposal.p_zero_freq=unifrnd(0,obj.prior.p_zero_freq_max);
                   case 2
                       found=false;
                       while ~found
                            temp=normrnd(proposal.p_zero_freq,obj.prior.p_zero_freq_max/obj.Lambda_2);
                            if temp>=0 && temp<=obj.prior.p_zero_freq_max
                                proposal.p_zero_freq=temp;
                                found=true;
                            end
                       end
               end
               MH_ratio=1;
               for j=1:length(proposal.freq_values)
                    if proposal.freq_values(j)>0
                        MH_ratio = MH_ratio*(1-proposal.p_zero_freq)/(1-old.p_zero_freq);
                    else
                        MH_ratio = MH_ratio*proposal.p_zero_freq/old.p_zero_freq;
                    end
               end
           else
               MH_ratio=1;
           end
        end
        function [proposal,MH_ratio]=propose_lambda_switches(obj,subtype,proposal)
            old=proposal;
            switch subtype
                case 1
                    proposal.lambda_switches=unifrnd(obj.prior.lambda_switches_min,obj.prior.lambda_switches_max);
                case 2
                    found=false;
                    while ~found
                        temp=normrnd(proposal.lambda_switches,(obj.prior.lambda_switches_max-obj.prior.lambda_switches_min)/obj.Lambda_2);
                        if temp>=obj.prior.lambda_switches_min && temp<=obj.prior.lambda_switches_max
                            proposal.lambda_switches=temp;
                            found=true;
                        end
                    end
            end
            MH_ratio = poisspdf(length(proposal.freq_switches),proposal.lambda_switches)/poisspdf(length(proposal.freq_switches),old.lambda_switches);
        end
        function [proposal,MH_ratio]=propose_freq_mean(obj,subtype,proposal)
            old=proposal;
            switch subtype
                case 1
                    proposal.freq_mean=unifrnd(obj.prior.freq_mean_min,obj.prior.freq_mean_max);
                case 2
                    found=false;
                    while ~found
                        temp=normrnd(proposal.freq_mean,(obj.prior.freq_mean_max-obj.prior.freq_mean_min)/obj.Lambda_2);
                        if temp>=obj.prior.freq_mean_min&& temp<=obj.prior.freq_mean_max
                            proposal.freq_mean=temp;
                            found=true;
                        end
                    end
            end
            MH_ratio = 1;
            for j=1:length(proposal.freq_values)
                if proposal.freq_values(j)>0
                MH_ratio = MH_ratio*invgampdf(proposal.freq_values(j),proposal.freq_mean,proposal.freq_alpha)/...
                    invgampdf(proposal.freq_values(j),old.freq_mean,old.freq_alpha);
                end
            end
        end
        function [proposal,MH_ratio]=propose_freq_alpha(obj,~,proposal)
            old=proposal;
            switch unidrnd(2)
                case 1
                    proposal.freq_alpha=unifrnd(obj.prior.freq_alpha_min,obj.prior.freq_alpha_max);
                case 2
                    found=false;
                    while ~found
                        temp=normrnd(proposal.freq_alpha,(obj.prior.freq_alpha_max-obj.prior.freq_alpha_min)/obj.Lambda_2);
                        if temp>=obj.prior.freq_alpha_min&& temp<=obj.prior.freq_alpha_max
                            proposal.freq_alpha=temp;
                            found=true;
                        end
                    end
            end
            MH_ratio = 1;
            for j=1:length(proposal.freq_values)
                if proposal.freq_values(j)>0
                MH_ratio = MH_ratio*invgampdf(proposal.freq_values(j),proposal.freq_mean,proposal.freq_alpha)/...
                    invgampdf(proposal.freq_values(j),old.freq_mean,old.freq_alpha);
                end
            end
        end
        function [proposal,MH_ratio]=propose_T_init(obj,subtype,proposal)
            old=proposal;
            switch unidrnd(2)
                case 1
                    proposal.T_init=unifrnd(obj.prior.T_init_min,obj.prior.T_init_max);
                case 2
                    found=false;
                    while ~found
                        temp=normrnd(proposal.T_init,(obj.prior.T_init_max-obj.prior.T_init_min)/obj.Lambda_2);
                        if temp>=obj.prior.T_init_min && temp<=obj.prior.T_init_max
                            proposal.T_init=temp;
                            found=true;
                        end
                    end
            end
            [proposal,MH_ratio]=obj.add_slip_history(old,proposal,subtype);
        end
        function [proposal,MH_ratio]=propose_parameters(obj,~,proposal)
            parameter_names=fieldnames(obj.parameters);
            i=unidrnd(length(obj.var_parameters));
            parameter_name=parameter_names{obj.var_parameters(i)};
            [proposal.(parameter_name),MH_ratio]=...
            obj.parameters.(parameter_name)(proposal.(parameter_name));
        end
        function [proposal,MH_ratio]=propose_S(obj,subtype,proposal)
           old=proposal;
           proposal.S=unifrnd(0,1);
           [proposal,MH_ratio]=obj.add_slip_history(old,proposal,subtype);
        end
        function [proposal,MH_ratio]=propose_bm(obj,subtype,proposal)
            old=proposal;
            interval_exponent=poissrnd(obj.lambda_interval_bm);
            interval_exponent=min(interval_exponent,log2(obj.N_bm-1)-2);
            c_intervals=2^(interval_exponent);
            i_interval=unidrnd(c_intervals);
            n_init=(i_interval-1)*(obj.N_bm-1)/c_intervals+1;
            n_end=i_interval*(obj.N_bm-1)/c_intervals+1;
            if i_interval<c_intervals   
                proposal.bm=obj.generate_bm(proposal.bm,n_init,n_end);
            else
                proposal.bm=obj.generate_bm(proposal.bm,n_init);  
            end
            [proposal,MH_ratio]=obj.add_slip_history(old,proposal,subtype);
        end
        function [proposal,MH_ratio]=propose_jump_sizes(obj,subtype,proposal)
            i_pre = max_such_that(proposal.times,@(x) x<proposal.T_init);%Last eroded earthquake
            after_jumps = proposal.jumps(i_pre+1:end); %Survived earthquakes
            N=length(after_jumps);
            if subtype==1%global
                MH_ratio=1;
                J_remain=obj.J_total;
                if (N*obj.d_min<J_remain)&&(N*obj.d_max>J_remain);
                    jumps=randfixedsum(N,1,J_remain,obj.d_min,obj.d_max); 
                    proposal.jumps=[unifrnd(obj.d_min,obj.d_max,i_pre,1);jumps];
                else
                    MH_ratio=0;
                end
            elseif subtype==2%local
                MH_ratio=1;
                N_redo=min(N-1,geornd(0.5)+1);
                i_redo=randsample(N,N_redo);
                sum_new_after_jumps=sum(after_jumps(i_redo));
                new_after_jumps=randfixedsum(N_redo,1,sum_new_after_jumps,obj.d_min,obj.d_max);
                after_jumps(i_redo)=new_after_jumps;
                proposal.jumps = [proposal.jumps(1:i_pre);after_jumps];
            end 
        end
        function [proposal,MH_ratio]=propose_freq_switches(obj,subtype,proposal)%Move one switch
            MH_ratio=1;
            old=proposal;
            if ~isempty(proposal.freq_switches)
                i=unidrnd(length(proposal.freq_switches));
                if i==1
                    T_below=obj.T_min;
                else
                    T_below=proposal.freq_switches(i-1);
                end
                if i==length(proposal.freq_switches)
                    T_above=obj.T_max;
                else
                    T_above=proposal.freq_switches(i+1);
                end
                T_old=proposal.freq_switches(i);
                T_new=unifrnd(T_below,T_above);
                proposal.freq_switches(i)=T_new;
                switch unidrnd(3)
                    case 1%keep same amount of earthquakes on each side of switch
                        L_1=T_new-T_below;
                        L_2=T_above-T_new;
                        L_1_old=T_old-T_below;
                        L_2_old=T_above-T_old;
                        proposal.freq_values(i)=L_1_old*old.freq_values(i)/L_1;
                        proposal.freq_values(i+1)=L_2_old*old.freq_values(i+1)/L_2;
                        MH_ratio = MH_ratio*invgampdf(proposal.freq_values(i),proposal.freq_mean,proposal.freq_alpha)/invgampdf(old.freq_values(i),old.freq_mean,old.freq_alpha);
                        MH_ratio = MH_ratio*invgampdf(proposal.freq_values(i+1),proposal.freq_mean,proposal.freq_alpha)/invgampdf(old.freq_values(i+1),old.freq_mean,old.freq_alpha);
                    case 2%redraw earthquake frequences
                        for j=1:2
                            if unifrnd(0,1)<proposal.p_zero_freq
                                proposal.freq_values(i-1+j)=0;
                            else
                                proposal.freq_values(i-1+j)=invgamrnd(proposal.freq_mean,proposal.freq_alpha);
                                %proposal.freq_values(i-1+j)=unifrnd(obj.prior.freq_min,proposal.freq_max);
                            end
                        end
                end 
                [proposal,MH_factor]=obj.add_slip_history(old,proposal,subtype);
                MH_ratio=MH_factor*MH_ratio;
                if ~obj.is_feasible_freq(proposal.freq_values,proposal)
                	MH_ratio=0;
                end
            else
                MH_ratio=0;
            end
        end
        function [proposal,MH_ratio]=propose_freq_values(obj,subtype,proposal)
            MH_ratio=1;
            old=proposal;
            N=length(proposal.freq_values);
            N_redo=min(N,geornd(0.5)+1);
            i_redo=randsample(N,N_redo);
            for i = i_redo'
                if unifrnd(0,1)<proposal.p_zero_freq
                    proposal.freq_values(i)=0;
                else
                    switch unidrnd(2)
                        case 1
                            proposal.freq_values(i) = invgamrnd(obj.state.freq_mean,obj.state.freq_alpha);
                            %proposal.freq_values(i)=unifrnd(obj.prior.freq_min,proposal.freq_max);%REP
                        case 2
                            found=false;
                            while ~found
                                temp=normrnd(proposal.freq_values(i),proposal.freq_mean/2/obj.Lambda_2);
                                if temp>=0
                                    found=true;
                                    proposal.freq_values(i)=temp;%REP ADD MH_RATIO
                                    MH_ratio = MH_ratio*invgampdf(proposal.freq_values(i),proposal.freq_mean,proposal.freq_alpha)/...
                                        invgampdf(old.freq_values(i),old.freq_mean,old.freq_alpha);
                                end
                            end
                    end
                end
            end 
            if ~obj.is_feasible_freq(proposal.freq_values,proposal)
                MH_ratio=0;
            end
            [proposal,MH_factor]=obj.add_slip_history(old,proposal,subtype);
            MH_ratio=MH_ratio*MH_factor;
        end
        function [proposal,MH_ratio]=propose_number_freq_switches(obj,subtype,proposal)
            if subtype==1%global
                %MH_ratio=1;
                N=poissrnd(proposal.lambda_switches);
                freq_switches_rand=unifrnd(obj.T_min,obj.T_max,[N,1]);
                proposal.freq_switches=sort(freq_switches_rand);
                proposal.freq_values=zeros(N+1,1);
                for j=1:N+1
                    if unifrnd(0,1)>=proposal.p_zero_freq
                        proposal.freq_values(j) = invgamrnd(obj.state.freq_mean,obj.state.freq_alpha);
                        %proposal.freq_values(j)=unifrnd(obj.prior.freq_min,proposal.freq_max);%REP
                    end
                end
                [proposal,MH_ratio]=obj.add_slip_history(proposal,proposal,subtype);
            else%local
                l=unifrnd(0,obj.T_total);
                T_below=unifrnd(obj.T_min-l,obj.T_max-l);
                change_first_freq=(T_below<obj.T_min);
                T_below=max(T_below,obj.T_min);
                T_above=T_below+l;
                N=poissrnd(l/obj.T_total*proposal.lambda_switches);
                new_switches=sort(unifrnd(T_below,T_above,[N,1]));
                new_values=zeros(N,1);
                for j=1:N
                    if unifrnd(0,1)>=proposal.p_zero_freq
                        new_values(j) = invgamrnd(obj.state.freq_mean,obj.state.freq_alpha);
                        %new_values(j)=unifrnd(obj.prior.freq_min,proposal.freq_max);%REP
                    end
                end
                i_below=max_such_that(proposal.freq_switches,@(x)x<T_below);
                i_above=min_such_that(proposal.freq_switches,@(x)x>T_above);
                if isempty(i_below)
                    i_below=0;
                end
                proposal.freq_switches=[proposal.freq_switches(1:i_below);new_switches;proposal.freq_switches(i_above:end)];
                proposal.freq_values=[proposal.freq_values(1:i_below+1);new_values;proposal.freq_values(i_above+1:end)];
                if change_first_freq
                    if unifrnd(0,1)>=proposal.p_zero_freq
                        proposal.freq_values(1)=invgamrnd(obj.state.freq_mean,obj.state.freq_alpha);
                    %proposal.freq_values(1)=unifrnd(obj.prior.freq_min,proposal.freq_max);%REP
                    else
                        proposal.freq_values(1)=0;
                    end
                end
                assert(length(proposal.freq_values)==(length(proposal.freq_switches)+1))
                [proposal,MH_ratio]=obj.add_slip_history(proposal,proposal,subtype);
            end
        end
        
        function [proposal,MH_ratio]=propose_freq_interlude(obj,subtype,proposal)
            old=proposal;
            N_freq_intervals=length(proposal.freq_switches)+1;
            N_freq_switches=N_freq_intervals-1;
            can_remove=(N_freq_intervals>2);
            addremove=unidrnd(1+can_remove);
            if addremove==1 %add switches
                interval=unidrnd(N_freq_intervals);
                if isempty(proposal.freq_switches)||interval==1
                    T_below=obj.T_min;
                else
                    T_below=proposal.freq_switches(interval-1);
                end
                if isempty(proposal.freq_switches)||interval==N_freq_intervals
                    T_above=obj.T_max;
                else
                    T_above=proposal.freq_switches(interval);
                end
                i=interval-1;
                [T_new,~]=sort(unifrnd(T_below,T_above,2,1));
                T_new_1=T_new(1);
                T_new_2=T_new(2);
                proposal.freq_switches=[proposal.freq_switches(1:interval-1);...
                                        T_new_1;...
                                        T_new_2;...
                                        proposal.freq_switches(interval:end)];
                MH_ratio=poisspdf(N_freq_switches+2,proposal.lambda_switches)/poisspdf(N_freq_switches,proposal.lambda_switches);%prior ratio
                new_freqs=zeros(3,1);
                for j=1:3
                    if unifrnd(0,1)<proposal.p_zero_freq
                        new_freqs(j)=0;
                    else
                        new_freqs(j) = invgamrnd(obj.state.freq_mean,obj.state.freq_alpha);
                        %new_freqs(j)=unifrnd(obj.prior.freq_min,proposal.freq_max);%rep
                    end
                end
                proposal.freq_values=[proposal.freq_values(1:i);...
                                        new_freqs;
                                        proposal.freq_values(i+2:end)];
                if ~obj.is_feasible_freq(proposal.freq_values,proposal)
                    MH_ratio=0;%prior ratio, part 2
                end
                if ~can_remove
                    MH_ratio=MH_ratio/2;% proposal ratio, part 2
                end
                [proposal,MH_factor]=obj.add_slip_history(old,proposal,subtype);
                MH_ratio=MH_ratio*MH_factor;%
            else %remove switch
                interval=unidrnd(N_freq_intervals-2)+1;
                if subtype==1
                    if unifrnd(0,1)<proposal.p_zero_freq
                        new_freq=0;
                    else
                        new_freq = invgamrnd(obj.state.freq_mean,obj.state.freq_alpha);
                        %new_freq=unifrnd(obj.prior.freq_min,proposal.freq_max);%REP
                    end
                    proposal.freq_values(interval-1)=new_freq;
                end%else, just keep the frequency that is currently before the interlude-to-be-removed
                proposal.freq_values(interval:interval+1)=[];
                proposal.freq_switches(interval-1:interval)=[];
                MH_ratio=poisspdf(N_freq_switches-2,proposal.lambda_switches)/poisspdf(N_freq_switches,proposal.lambda_switches);%prior ratio
                if ~obj.is_feasible_freq(proposal.freq_values,proposal)
                    MH_ratio=0;%prior ratio, part 2. necessary?
                end
                if N_freq_intervals==2
                    MH_ratio=MH_ratio*2;% proposal ratio, part 2
                end
                [proposal,MH_factor]=obj.add_slip_history(old,proposal,subtype);
                MH_ratio=MH_ratio*MH_factor;%
            end
        end
        
        function [proposal,MH_ratio]=propose_tau(obj,subtype,proposal)
            old=proposal;
            delta_tau=obj.prior.tau_max-obj.prior.tau_min;
            switch unidrnd(2)
                case 1
                    proposal.tau=unifrnd(obj.prior.tau_min,obj.prior.tau_max);
                case 2
                    found = false;
                    while ~found
                        temp=normrnd(proposal.tau,delta_tau/obj.Lambda_2);
                        if temp<=obj.prior.tau_max && temp>=obj.prior.tau_min
                            found=true;
                            proposal.tau=temp;
                        end
                    end
            end
            if obj.is_feasible_tau(proposal.tau)
                MH_ratio=1;
            else
                MH_ratio=0;
            end
            [proposal,MH_factor]=obj.add_slip_history(old,proposal,subtype);
            MH_ratio=MH_ratio*MH_factor;
        end
        function bm=generate_bm(obj,bm,n_init,n_end)
            if nargin<4 || isempty(n_end) || n_end==obj.N_bm%Brownian motion
                temp=normrnd(0,sqrt(1/obj.N_bm),obj.N_bm-n_init,1);
                bm=[bm(1:n_init);bm(n_init)+cumsum(temp)];
            else %Brownian bridge
               if (n_end-n_init)>=2
                   bm_init=bm(n_init);
                   bm_end=bm(n_end);
                   n_mean=floor((n_init+n_end)/2);
                   bm_mean=normrnd(bm_init+(bm_end-bm_init)*(n_mean-n_init)/(n_end-n_init),sqrt((n_end-n_mean)*(n_mean-n_init)/(n_end-n_init)*(1/obj.N_bm)));
                   bm(n_mean)=bm_mean;
                   bm=obj.generate_bm(bm,n_init,n_mean);
                   bm=obj.generate_bm(bm,n_mean,n_end);
               end
            end
        end
        function [event_times,time_after,X]=find_times(obj,state)
            if obj.is_feasible_freq(state.freq_values,state)
                X=zeros(obj.N_bm,1);
                T=obj.T_min+((1:obj.N_bm)'-1)*obj.T_total/(obj.N_bm-1);
                X(1)=state.S;
                event_times=[];
                i_freq=1;
                n_switches=length(state.freq_switches);
                for n=2:obj.N_bm
                    if n_switches>=i_freq && T(n)>state.freq_switches(i_freq)
                        i_freq=i_freq+1;
                    end
                    freq=state.freq_values(i_freq);
                    bm_drift=freq;
                    bm_diffusivity=sqrt(bm_drift)*state.tau;
                    dt=obj.T_total/obj.N_bm;
                    X(n)=X(n-1)+bm_drift*dt+bm_diffusivity*(state.bm(n)-state.bm(n-1))*sqrt(obj.T_total);
                    if X(n)<0
                    elseif X(n)>=1
                       event_times=[event_times;T(n)*ones(floor(X(n)),1)];
                       X(n)=X(n)-floor(X(n));
                    end
                end
                time_after=0;
                Xafter=X(obj.N_bm);
                freq=state.freq_values(end);
                found=false;
                dt=10;
                while ~found
                    next_switch=exprnd(1/state.lambda_switches*obj.T_total);
                    N=next_switch/dt;
                    rd=normrnd(0,sqrt(dt),ceil(N),1);
                    for n=1:N
                        time_after=time_after+dt;
                        bm_drift=freq;
                        bm_diffusivity=sqrt(bm_drift)*state.tau;
                        Xafter=Xafter+bm_drift*dt+bm_diffusivity*rd(n);
                        if Xafter>=1||time_after>obj.cutoff
                           found=true;
                           break
                        end
                    end
                    if unifrnd(0,1)<state.p_zero_freq
                        freq=0;
                    else
                        freq=invgamrnd(obj.state.freq_mean,obj.state.freq_alpha);   
                    end
                end
            else
                event_times=-Inf;%will be discarded
                X=-Inf;%will be discarded
                time_after=-Inf;
            end
        end
        function B=is_feasible_N(~,N)
            B=(N>0);
        end
        function B=is_feasible_freq(~,d,~)
            B = all(d>=0);
        end
        function B=is_feasible_tau(obj,tau)
           B = (tau<=obj.prior.tau_max)&&(tau>=obj.prior.tau_min); 
        end
        function B=is_feasible_slip_history(obj,proposal)
            i_pre = max_such_that(proposal.times,@(x) x<proposal.T_init);%Last eroded earthquake
            pre_jumps = proposal.jumps(1:i_pre);
            after_jumps = proposal.jumps(i_pre+1:end);
            B=all(pre_jumps<=obj.d_max)&&all(pre_jumps>=obj.d_min)&&...
                all(after_jumps<=obj.d_max)&&all(after_jumps>=obj.d_min)&&...
                abs(sum(after_jumps)-obj.J_total)<1e-4;
        end
        function [proposal,MH_factor]=add_slip_history(obj,old,proposal,subtype)
            [proposal.times,proposal.time_after]=obj.find_times(proposal);
            new_i_pre = max_such_that(proposal.times,@(x) x<proposal.T_init);%Last eroded earthquake
            new_times = proposal.times(new_i_pre+1:end);%Only those that survived
            N_new=length(new_times);
            if N_new>0
                MH_factor=1;
                if subtype==1%Global
                    J_remain=obj.J_total;
                    if (N_new*obj.d_min<J_remain)&&(N_new*obj.d_max>J_remain)
                        jumps=randfixedsum(N_new,1,J_remain,obj.d_min,obj.d_max);
                    else
                        MH_factor = 0;
                        jumps=-Inf;%will be discarded
                    end
                elseif subtype==2%local   
                    old_i_pre = max_such_that(old.times,@(x) x<proposal.T_init);%Last eroded earthquake
                    old_times = old.times(old_i_pre+1:end); %Only those that survived
                    old_jumps = old.jumps(old_i_pre+1:end);
                    N_old=length(old_times);
                    if N_new==N_old
                        jumps=old_jumps;
                    elseif N_new>N_old
                        J_remain=obj.J_total;
                        sum_new_jumps=(N_new-N_old)*obj.d_min+(J_remain-(N_new-N_old)*obj.d_min)*betarnd(N_new-N_old,N_old);
                        if ((N_new-N_old)*obj.d_min<sum_new_jumps)&&((N_new-N_old)*obj.d_max>sum_new_jumps)
                            new_jumps=randfixedsum(N_new-N_old,1,sum_new_jumps,obj.d_min,obj.d_max);
                        else%Could happen with too many new jumps
                            new_jumps=-Inf*ones(N_new-N_old,1);%will be discarded
                        end
                        sum_old_jumps=J_remain-sum_new_jumps;
                        multiply_old=sum_old_jumps/sum(old_jumps);
                        temp=multiply_old*old_jumps;
                        i_insert=unidrnd(N_old+1);%Insert before this index
                        jumps=[temp(1:i_insert-1);new_jumps;temp(i_insert:N_old)];%
                    elseif N_new<N_old
                        J_remain=obj.J_total;
                        i=unidrnd(N_new+1);%starting here remove N_old-N_new
                        temp=[old_jumps(1:i-1);old_jumps(i+N_old-N_new:end)];
                        jumps=J_remain/sum(temp)*temp;
                    end                    
                end
            else
                jumps=zeros(0,1);
                MH_factor=0;
            end              
            proposal.jumps=[unifrnd(obj.d_min,obj.d_max,new_i_pre,1);jumps];
            if ~obj.is_feasible_slip_history(proposal)
                MH_factor=0;
            end
        end
        function plot_tension(obj,state)
            [times,X]=obj.find_times(state);
            hold off
            plot(linspace(obj.T_min,obj.T_max,length(X)),X)
            hold on
            scatter(times',zeros(1,length(times)),max(state.jumps,1),'filled')
            hold on 
            plot(state.freq_switches,zeros(1,length(state.freq_switches)),'bo');
        end
    end
end
