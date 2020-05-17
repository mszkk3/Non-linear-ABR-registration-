function [lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx,varargin)
%xtractlatamp Extraction of ABR wave latencies and amplitudes
%   
%   Syntax:
%       [lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx)
%       [lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx,lag)
%       [lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx,lag/[],tgt)
%       [lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx,lag/[],[],TGT)
%
%   Deacription:
%       [lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx) extracts latencies
%       and amplitudes of an arbitrary number of waves using the warped ABRs,
%       wabr, and associated time warping functions, tstar, generated by 
%       nlcurvereg.m. Structural latencies are picked in the structual 
%       grand-average response and the pre-alignment lags are assumed 
%       to be zero. 
%       wabr ==> are the warped ABRs created by nlcurvereg.m.
%       tstar ==> are the corresponding warping functions.
%       tn ==> is the normalized time axis created by nlcurvereg.m.
%       time ==> is the original time sample array also used in preproc.m
%       and nlcurvereg.m.
%       tidx ==> is the index of time samples within the warping time range
%       created by preproc.m.
%       lat ==> array or cell array of individual peak and trough latencies 
%       for each wave pcicked. lat or each of its cells will have dimensions
%       N (number of subjects) by 2 (peak/trough) by M, where M is the number 
%       of waves picked.
%       amp ==> array or cell array of individual peak and trough
%       amplitudes with the same dimensions as lat. 
%       gixd ==> indices of structural peak and trough latencies within
%       tim. 
%
%       [lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx,lag) same as above,
%       except that pre-alignment lags, lag - generated by preproc.m are
%       taken into account when deriving the wave latencies. 
%
%       [lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx,lag/[],tgt) same above,
%       except that structural peak and trough latencies are picked in a
%       predefined structural target response, tgt. If lag = [],
%       pre-alignment lags are assumed to be zero. 
%
%       [lat,amp,gidx] = xtractlatamp(wabr,tstar,tim,tidx,lag/[],[],TGT) if abr/wabr is a 
%       cell array with C cells and TGT is a scalar integer between 1 and
%       C, the structural latencies will be picked in the abverage response for the 
%       TGTth condition (i.e., mean(wabr{TGT})). If lag = [], the pre-alignment 
%       lags are assumed to be zero. 

   %% Read input.
    vars = {'lag' 'tgt' 'TGT'};
    for I = 1:nargin-4
        eval(sprintf('%s = varargin{%d};',vars{I},I))
    end

    %% Create normalized time axis.
    T = numel(tim);    
    t = tim(tidx);
    tn = (t-t(1))/(t(end)-t(1)); 

    %% Pick wave peaks and troughs. 
    if iscell(wabr)
        C = numel(wabr);
        N = size(wabr{1},1);
    else
        C = 1;
        N = size(wabr,1);
    end
    
    if ~exist('lag','var')||isempty(lag)
        if C>1
            lag = zeros(1,C);
        else
            lag = zeros(1,N);
        end
    end

    if ~exist('tgt','var')||isempty(tgt)
        if ~exist('TGT','var')||isempty(TGT)
            if iscell(wabr)
                tmp = cell2mat(wabr'); 
            else
                tmp = wabr;
            end
        else
            tmp = wabr{TGT};
        end
        sa = mean(tmp(:,tidx)); sa = sa-mean(sa);
    else
        if numel(tgt)==T, sa = tgt(tidx);  end
        sa = sa-mean(sa);        
    end
    dsa = sign(diff([sa sa(end)])); dsa = diff([dsa(1) dsa]);             
        
    H = figure; set(gcf,'Name','Pick waves'), hold on
    set(gca,'XLim',[t(1) t(end)])
    xlabel('Time (ms)'), ylabel('Amplitude (\muV)')
    line(xlim,zeros(1,2),'Color','k')
    plot(t,sa,'k-','LineWidth',2)
    pidx = find(dsa<0); nidx = find(dsa>0);
    for I = 1:numel(pidx)
        plot(t(pidx(I)),sa(pidx(I)),'r+')
    end
    for I = 1:numel(nidx)
        plot(t(nidx(I)),sa(nidx(I)),'b+')
    end       
    
    gidx = [];
    WavNum = 1;
    while true
        fprintf(1,'\nPick peak of wave #%d [or press <return> to exit]\n',WavNum);
        figure(H), xy = ginput(1); 
        if ~isempty(xy)
            [~,IDX1] = min(abs(t-xy(1)));
            L1 = line(repmat(t(IDX1),1,2),ylim,'Color','m','LineStyle','--');
            [~,IDX2] = min(abs(t(pidx)-t(IDX1)));
            L2 = line(repmat(t(pidx(IDX2)),1,2),ylim,'Color','r','LineStyle','--');
            Move = input('Move to nearest peak ([y]/n)?','s');
            if strcmp(Move,'n')
                delete(L2)
                PIDX = IDX1;
                line(repmat(t(PIDX),1,2),ylim,'Color','r','LineStyle','--');
            else
                delete(L1)
                PIDX = pidx(IDX2);
            end

            fprintf(1,'Pick trough of wave #%d [or press <return> to excit]\n',WavNum);
            figure(H), xy = ginput(1); 
            if ~isempty(xy)                
                [~,IDX1] = min(abs(t-xy(1)));
                L1 = line(repmat(t(IDX1),1,2),ylim,'Color','c','LineStyle','--');
                [~,IDX2] = min(abs(t(nidx)-t(IDX1)));
                L2 = line(repmat(t(nidx(IDX2)),1,2),ylim,'Color','b','LineStyle','--');
                Move = input('Move to nearest trough ([y]/n)?','s');
                if strcmp(Move,'n')
                    delete(L2)
                    TIDX = IDX1;
                    line(repmat(t(TIDX),1,2),ylim,'Color','b','LineStyle','--');
                else
                    delete(L1)
                    TIDX = nidx(IDX2);
                end
                gidx = [gidx;[PIDX TIDX]];
                WavNum = WavNum+1;
            else
                break
            end
        else
            break
        end
    end
    delete(H)
    
    %% Derive individual latencies and amplitudes.
    if ~isempty(gidx)
        M = size(gidx,1);
        if iscell(wabr)
            lat = cell(size(wabr));
            amp = cell(size(wabr));
            for I = 1:C
                lat{I} = zeros(N,2,M);
                amp{I} = zeros(N,2,M);
                for II = 1:M
                    lat{I}(:,:,II) = interp1(tn,tstar{I}',tn(gidx(II,:)))'*(tim(end)-tim(1))+tim(1);
                    amp{I}(:,:,II) = wabr{I}(:,tidx(gidx(II,:)));
                end
            end
        else
            lat = zeros(N,2,M);
            amp = zeros(N,2,M);
            for I = 1:M
                lat(:,:,I) = interp1(tn,tstar',tn(gidx(I,:)))'*(tim(end)-tim(1))+tim(1);
                amp(:,:,I) = wabr(:,tidx(gidx(I,:)));
            end
        end
        
        gidx = tidx(gidx);
    end
end
            
            
