%%
%% Automatic axon extraction
% USAGE
%
%   [branch] = axon_velocity_neurons_auto_5(neurons,ntk,param,plot_on)
%   neurons     sorted neuron template struct including
%               struct with fields:
%                 templates_2   %templates (length*no.channels)
%                  spike_time   %spike time
%                 channel_all   %init channel no.
%                    valid_el   %valid electrodes after pre-cleaning
%   plot_on     for plotting
%   ntk     struct with array info
%   param   extraction parameters
%             param.min_points=10;% minimum number of points for each branch
%             param.jump_space=80;% searching radius
%             param.jump_frame = 5;% no.frames to jump over
%             param.init_delay=5;% no.frames after init point
% OUTPUT
%   branch    %extracted branches
%             selected_axon %selected els
%             velocity, offset, error %fitting results
%             space,time %data for fitting

 function [branch] = axon_velocity_auto(neuron,ntk,param,plot_on)

%%
waveforms = neuron.templates_2(:,neuron.valid_el);% use only valid electrodes after pre-cleaning
x= ntk.y2(neuron.valid_el);
y=ntk.x2(neuron.valid_el);
%%
[~, peak_time] = min(waveforms);
init_ch=find(neuron.valid_el==neuron.channel_all);
init_frame = peak_time(init_ch);

%% Find last frames as starting points
last_axon_frames=find((peak_time>=init_frame+param.init_delay));
[~,ind_frame] = sort(peak_time(last_axon_frames),'descend');
   
%%
clear space
branch=cell(0,0);
selected_axon_all=[];
max_distance=[];
frame_last=[];
%% Main extraction part
for branch_no=1:length(last_axon_frames)
    last_ch = last_axon_frames(ind_frame(branch_no));% channel no for last frame
    branch(branch_no).selected_axon=[];
    clear space
    axon_ch=last_ch;
    axon_frame=peak_time(last_ch);% frame no for last frame
    selected_axon = [];
    count=0;
    while ~isempty(axon_ch)
        count = count+1;
        selected_axon(count) = axon_ch;
        next_axon_ch=[];
        for j=1:param.jump_frame
            clear space
            next_frame_ch = find(peak_time==axon_frame-j);% find previous frame
            if ~isempty(next_frame_ch)&((axon_frame-j)>peak_time(init_ch)+param.init_delay)
                for k = 1:length(next_frame_ch)
                    space(k) = sqrt((x(axon_ch)-x(next_frame_ch(k)))^2 + (y(axon_ch)-y(next_frame_ch(k)))^2);% calculate distance
                end
                next_frame_ch = next_frame_ch(find(space<param.jump_space*j));% find el within radius
                if ~isempty(next_frame_ch)
                    next_amp = min(waveforms(:,next_frame_ch));
                    [~,ind]=sort(next_amp,'ascend');
                    next_axon_ch = next_frame_ch(ind(1));% choose el with largest amplitude
                    break;
                end
            end
        end
        axon_ch=next_axon_ch;% new start point for seaching
        axon_frame = peak_time(axon_ch);
    end
    selected_axon=flip((selected_axon));% reverse the selected els
    %% integrate distances for all branches
    for i =1:size(selected_axon,2)
        selected(i) = selected_axon(i);
        if i==1
                space(i) = 0;
        else
                diff_space(i) = sqrt((x(selected(i))-x(selected(i-1)))^2 + (y(selected(i))-y(selected(i-1)))^2);
                space(i) = space(i-1)+diff_space(i);
        end
    end
    branch(branch_no).selected_axon=selected_axon;
    branch(branch_no).time = (peak_time(selected_axon)-peak_time(selected_axon(1)))/ntk.sr;
    selected_axon_all = [selected_axon_all selected_axon];
    max_distance(branch_no) = max(space);
    frame_last(branch_no)=last_axon_frames(ind_frame(branch_no));   
    clear space
end
clear space

%% Sort and clean branches
branch_2 = cell(0,0);
ind = [1:1:length(branch)];
velocity=[];
error=[];
count=1;
selected_axon_2=[];
for branch_no = 1:length(ind)
    % check if the selected el belongs to another branch
    for j = 1:length(branch(ind(branch_no)).selected_axon)
        if ~isempty(find(selected_axon_2==branch(ind(branch_no)).selected_axon(j)))
           branch(ind(branch_no)).selected_axon(j)=NaN;%remove existing els
        end
    end
    branch(ind(branch_no)).selected_axon = branch(ind(branch_no)).selected_axon(~isnan(branch(ind(branch_no)).selected_axon));
    %% output only branches with more than x points
    if length(branch(ind(branch_no)).selected_axon)>= param.min_points
        selected_axon = branch(ind(branch_no)).selected_axon;
         for i =1:size(selected_axon,2)
            selected(i) = selected_axon(i);
            if i==1
                    space(i) = 0;
            else
                    diff_space(i) = sqrt((x(selected(i))-x(selected(i-1)))^2 + (y(selected(i))-y(selected(i-1)))^2);
                    space(i) = space(i-1)+diff_space(i);
            end
        end
        [fit1,gof]=fit((peak_time(selected_axon)'-peak_time(selected_axon(1)))/ntk.sr*1e3, space','poly1');
        coeffs(1) = fit1.p1;
        coeffs(2) = fit1.p2;        
        branch_2{count}.selected_axon=selected_axon;
        branch_2{count}.velocity = coeffs(1);
        branch_2{count}.offset = coeffs(2);
        branch_2{count}.error=gof.adjrsquare;
        branch_2{count}.space = space;
        branch_2{count}.time = (peak_time(selected_axon)-peak_time(selected_axon(1)))/ntk.sr*1e3;
        selected_axon_2 = [selected_axon_2 (branch_2{count}.selected_axon)];% register all selected els for cleaning later
        velocity(count)=branch_2{count}.velocity;
        error(count) = branch_2{count}.error;
        count= count+1;
        clear space selected_axon selected
    end
end
branch = branch_2;
    

%% For plotting
if plot_on
[amp_aps delay_aps] = min(neuron.templates_2);
for i=1:1403
        for ii=1:722
            map_amp(i,ii) = amp_aps(ntk.channel_no2(i,ii));
            if find(ntk.channel_no2(i,ii)==neuron.valid_el)
                map_delay(i,ii) = (delay_aps(ntk.channel_no2(i,ii))-init_frame)/ntk.sr*1e3;
            else
                map_delay(i,ii)=Inf;
            end
        end
    end
end

[val ind]=sort(velocity);
branch = branch(ind);
figure('color','w','position',[100 200 1600 280]);
%% amplitude map
subplot(1,4,1)
imagesc(ntk.yona,ntk.xona,map_amp')
axis equal;
ylim([0 1786]) 
xlim([0 3450])
xticks([]);
yticks([]);
set (gca,'Xdir','reverse')
set (gca,'Ydir','default')
caxis([-20 0])
myColorMap=colormap(gca,flipud(parula(256)));
colormap(gca,myColorMap);
cb = colorbar; 
set(cb,'position',[.28 .32 .005 .15])
cb.AxisLocation='in';
cb.FontSize=8;
cb.Label.String = '\muV';

%% amp map with selected els
subplot(1,4,3)
imagesc(ntk.yona,ntk.xona,map_amp')
axis equal;
ylim([0 1786]) 
xlim([0 3450])
xticks([]);
yticks([]);
set (gca,'Xdir','reverse')
set (gca,'Ydir','default')
caxis([-20 0])
myColorMap=colormap(gca,flipud(parula(256)));
colormap(gca,myColorMap);
cb = colorbar; 
set(cb,'position',[.69 .32 .005 .15])
cb.AxisLocation='in';
cb.FontSize=8;
cb.Label.String = '\muV';
hold on
c = jet(size(branch,2));
for i = 1:size(branch,2)
    x_selected = x(branch{i}.selected_axon);
    y_selected = y(branch{i}.selected_axon);
    if i==1
        scatter(x_selected,y_selected,60,'.','MarkerEdgeColor',c(i,:))
    elseif i==length(branch)
        scatter(x_selected,y_selected,60,'.','MarkerEdgeColor',c(i,:))
    else
        scatter(x_selected,y_selected,60,'.','MarkerEdgeColor',c(i,:))
    end   
    hold on
end
axis equal;
ylim([0 1786]) 
xlim([0 3450])
ylim([0 1786]) 
xlim([0 3450])
xticks([]);
yticks([]);
set (gca,'Xdir','reverse')
set (gca,'Ydir','default')

%% delay map
subplot(1,4,2)
imagesc(ntk.yona,ntk.xona,(map_delay'))
hold on
axis equal;
ylim([0 1786]) 
xlim([0 3450])
xticks([]);
yticks([]);
set (gca,'Xdir','reverse')
set (gca,'Ydir','default')
cb = colorbar; 
set(cb,'position',[.485 .32 .005 .15])
cb.AxisLocation='in';
cb.FontSize=8;
cb.Label.String = 'ms';
myColorMap=colormap(gca,flipud(parula(256)));
myColorMap(255:256,:) = 1;
colormap(gca,myColorMap);
ylim([0 1786]) 
xlim([0 3450])
xticks([]);
yticks([]);
%% fitting
subplot('position',[0.75,0.3,0.1,0.4])
for i = 1:size(branch,2)  
    if i==1
        plot((peak_time(branch{i}.selected_axon)-peak_time(branch{i}.selected_axon(1)))/11.6,branch{i}.space/1000,'.','markersize',20,'color',c(i,:));
    elseif i ==length(branch)
        plot((peak_time(branch{i}.selected_axon)-peak_time(branch{i}.selected_axon(1)))/11.6,branch{i}.space/1000,'.','markersize',20,'color',c(i,:));
    else
        plot((peak_time(branch{i}.selected_axon)-peak_time(branch{i}.selected_axon(1)))/11.6,branch{i}.space/1000,'.','markersize',20,'color',c(i,:));
    end
    hold on
end

for i = 1:size(branch,2)
    xx = 0:.1:5;
    yy = (branch{i}.velocity*xx+branch{i}.offset)/1000;
    hold on;
    plot(xx,yy,'-','color',c((i),:))
    hold on
    xlim([0 5]);%max(branch{i}.space)])
    ylim([0 Inf]);
end
xlabel('Latency [ms]')
ylabel(['Distance[mm]'])
box off
text(0.5,3,['v = ',num2str(ceil(mean(velocity))),' \pm ',num2str(ceil(std(velocity))),'mm/s']);
text(0.5,2.5,['R-square = ',sprintf('%.3f',mean(error)),' \pm ',sprintf('%.3f',std(error))]);


end
