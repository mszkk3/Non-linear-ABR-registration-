load('chirpdata.mat')
who

%% Example 1: Warps the broadband chirp-evoked  ABRs (single condition 
% and 22 subjects) to to the subject-average response. Individual responses
% are not pre-aligned.

[obj,tidx] = preproc(abr{TGT},tim,wrng);
[wabr,tstar] = nlcurvereg(obj,tim,tidx,abr{TGT});
[lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx);

figure('Name','ABR registration - example 1')
subplot(1,2,1), hold on
t = tim(tidx);
set(gca,'XLim',[min(t) max(t)],'Ylim',[-0.6 0.6])
title('Original')
for I = 1:size(abr{TGT},1)
    plot(t,abr{TGT}(I,tidx),'Color',0.5*ones(1,3),'LineWidth',1)
end
plot(t,mean(abr{TGT}(:,tidx)),'k-','LineWidth',3)

subplot(1,2,2), hold on
set(gca,'XLim',[min(t) max(t)],'Ylim',[-0.6 0.6])
title('Aligned')
for I = 1:size(wabr,1)
    plot(t,wabr(I,tidx),'Color',0.5*[1 0 0]+0.5*ones(1,3),'LineWidth',1)
end
plot(t,mean(wabr(:,tidx)),'r-','LineWidth',3)

xlabel('Time (ms)')
ylabel('Amplitude (\muV')

%% Example 2: Warps all chirp-evoked ABRs (six high-pass-masking conditions
% and 22 subjects) to to the subject-average response for the broadband 
% condition. The subject-average responses are pre-aligned across conditions.

[obj,tidx,lag] = preproc(abr,tim,wrng,[],TGT,lag0,srng);
[wabr,tstar] = nlcurvereg(obj,tim,tidx,abr,lag,[],TGT);
[lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx,lag,[],TGT);

figure('Name','ABR registration - example 2')
t = tim(tidx);
for I = 1:numel(abr)
    subplot(2,3,I), hold on
    set(gca,'XLim',[min(t) max(t)],'Ylim',[-0.6 0.6])
    if I==2, title('Original'), end

    for II = 1:size(abr{I},1)
        plot(t,abr{I}(II,tidx),'Color',0.5*ones(1,3),'LineWidth',1)
    end
    plot(t,mean(abr{I}(:,tidx)),'k-','LineWidth',3)

    if I==5, xlabel('Time (ms)'), end
    if I==1||I==2, ylabel('Amplitude (\muV'), end
end

figure('Name','ABR registration - example 2')
for I = 1:numel(wabr)
    subplot(2,3,I), hold on
    set(gca,'XLim',[min(t) max(t)],'Ylim',[-0.6 0.6])
    if I==2, title('Aligned'), end

    for II = 1:size(wabr{I},1)
        plot(t,wabr{I}(II,tidx),'Color',0.5*[1 0 0]+0.5*ones(1,3),'LineWidth',1)
    end
    plot(t,mean(wabr{I}(:,tidx)),'r-','LineWidth',3)

    if I==5, xlabel('Time (ms)'), end
    if I==1||I==2, ylabel('Amplitude (\muV'), end
end

%% Example 3: Warps all chirp-evoked ABRs (six high-pass-masking conditions
% and 22 subjects) to to predefined structural target response (tgt). The 
% subject-average response are pre-aligned with tgt.

[obj,tidx,lag] = preproc(abr,tim,wrng,tgt,[],[],srng);
[wabr,tstar] = nlcurvereg(obj,tim,tidx,abr,lag,tgt);
[lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx,lag,tgt);

figure('Name','ABR registration - example 3')
t = tim(tidx);
for I = 1:numel(abr)
    subplot(2,3,I), hold on
    set(gca,'XLim',[min(t) max(t)],'Ylim',[-0.6 0.6])
    if I==2, title('Original'), end

    for II = 1:size(abr{I},1)
        plot(t,abr{I}(II,tidx),'Color',0.5*ones(1,3),'LineWidth',1)
    end
    plot(t,mean(abr{I}(:,tidx)),'k-','LineWidth',3)
    if I==TGT, plot(t,tgt(:,tidx),'b-','LineWidth',3), end

    if I==5, xlabel('Time (ms)'), end
    if I==1||I==2, ylabel('Amplitude (\muV'), end
end

figure('Name','ABR registration - example 3')
for I = 1:numel(wabr)
    subplot(2,3,I), hold on
    set(gca,'XLim',[min(t) max(t)],'Ylim',[-0.6 0.6])
    if I==2, title('Aligned'), end

    for II = 1:size(wabr{I},1)
        plot(t,wabr{I}(II,tidx),'Color',0.5*[1 0 0]+0.5*ones(1,3),'LineWidth',1)
    end
    plot(t,mean(wabr{I}(:,tidx)),'r-','LineWidth',3)
    if I==TGT, plot(t,tgt(:,tidx),'b-','LineWidth',3), end

    if I==5, xlabel('Time (ms)'), end
    if I==1||I==2, ylabel('Amplitude (\muV'), end
end
