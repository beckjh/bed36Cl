function VisualizeBED(input)
    % INPUT Name of Case Study or structure `results` from results.mat file
    % This file contains mostly GUI code. 
    % To change plots, look at BED_plots.m
    main_directory=fileparts(mfilename('fullpath'));
    cd(main_directory)
    addpath(genpath('.'))
    if ischar(input)
        input_file = ['CaseStudies/',input,'/results.mat'];    
        fprintf('Loading results from %s ...',input_file);
        tmp = load(input_file);
        results = tmp.results;
        fprintf(' done\n');
    else
        if length(fieldnames(input)) == 1
            results = input.results;
        else
            results = input;
        end
    end

    try
        close all
    catch
        figure(1)
        set(f, 'CloseRequestFcn',@(~,~)closereq())
        close all
    end
    f=figure(1);
    set(gcf,'Name','VisualizeBED','NumberTitle','off')
    set(f, 'Visible', 'off');
    set(f, 'menubar', 'none');
    set(f, 'pos', [0 0 235 150]);
    set(f, 'CloseRequestFcn',@close_request);
    plt =  uicontrol('Style', 'pushbutton', 'String', 'Plot',...
        'Position', [5 4 225 25]);
    uicontrol('Style','text','Position',[0 33 90 20],'String','Burn in (%)');
    brn = uicontrol('Parent',f,'Style','edit','Position',[116,33,50,23],...
        'String',10,'Keypressfcn',@key_press);
    t2 = uicontrol('Parent',f,'Style','edit','Position',[177,60,50,23],...
        'String',num2str(floor(results.settings.T_max/1000)),...
        'Keypressfcn',@key_press);
    t1 = uicontrol('Parent',f,'Style','edit','Position',[116,60,50,23],...
        'String',num2str(floor(results.settings.T_min/1000)),...
        'Keypressfcn',@key_press);
    uicontrol('Parent',f,'Style','text','Position',[167,70,10,10],...
        'String','--');
    uicontrol('Style','text','Position',[0 60 75 20],'String','Plot time'); 
	cnf = uicontrol('Parent',f,'Style','edit','Position',[116,90,50,23],...
        'String',90,...
        'Keypressfcn',@key_press);
    uicontrol('Style','text','Position',[2 90 110 20],'String','Confidence (%)');
    res = uicontrol('Parent',f,'Style','edit','Position',[116,120,50,23],...
        'String',1,...
        'Keypressfcn',@key_press);
    uicontrol('Style','text','Position',[5 120 110 20],'String','Resolution factor');
	distributed=0;
    T=@(es,ed) do_something(brn.String,'1',t1.String,t2.String,cnf.String,res.String);
    plt.Callback=T;
    f.Children=f.Children(end:-1:1);
    T(1,1);
    function key_press(~, eventdata, ~)
        if strcmp(eventdata.Key,'return')
           pause(0.01)
           T(1,1); 
        end
    end
    function close_request(~,~)
        c=allchild(0);
        for i = 1:length(c)
            delete(c(i));
        end
    end  
    function do_something(p_burn_in,level,plot_time_start,plot_time_end,confidence,resolution_factor)
        confidence = str2double(confidence);
        p_burn_in = str2double(p_burn_in)/100;
        plot_time = [1e3*str2double(plot_time_start),1e3*str2double(plot_time_end)];
        level = round(str2double(level));
        resolution_factor = str2double(resolution_factor);
        range = 1+results.settings.group_size*(level-1):results.settings.group_size*level;
        BED_plots(results,p_burn_in,plot_time,confidence,range,resolution_factor);
        if ~distributed && ishandle(2)
            distributed = 1;
            set(f,'Visible','off');
            distFig();
        end
        set(f,'Visible','on');
    end
end
