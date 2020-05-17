function [wabr,tstar] = nlcurvereg(obj,tim,tidx,abr,varargin)
%nlcurvereg non-linear curve registration of ABRs
%
%   Syntax:
%       [wabr,tstar] = nlcurvereg(obj,tim,tidx,abr)
%       [wabr,tstar] = nlcurvereg(obj,tim,tidx,abr,lag)
%       [wabr,tstar] = nlcurvereg(obj,tim,tidx,abr,lag/[],tgt)
%       [wabr,tstar] = nlcurvereg(obj,tim,tidx,abr,lag/[],[],TGT)
%
%   Description:
%       [wabr,tstar] = nlcurvereg(obj,tim,tidx,abr) non-linearly
%       alings the ABRs contained in abr using the average-target
%       registration procedure with two iterations, the PSDD fitting criterion 
%       (penalised sum of squared differences of response derivatives), and no 
%       warping roughness penalty (lambda = 0). The warping target is the 
%       grand-average response across all subjects and conditions (if using 
%       multiple conditions). It is assumed that the responses were not subjected 
%       to any pre-alignment (i.e., the pre-alignment lags were
%       zero).
%       abr ==> is either an array with dimensions N (number of subjects) by T
%       (number of time points), or a cell array with C cells, where C is
%       the number of conditions, and each cell has dimensions N by T.
%       tim ==> is an array of the T time measurement points in milliseconds. 
%       obj are the pre-processed ABRs, generated using preproc.m.
%       tidx ==> is the index of W time samples within the warping time range 
%       created by preproc.m.
%       wabr ==> is an N-by-T array or 1-by-C cell array (if abr is a cell array)  
%       containing the warped ABRs.
%       tstar ==> is an N-by-W array or 1-by-C cell array (if abr is a cell array)   
%       containing the time warping functions. 
%
%       [wabr,tstar] = nlcurvereg(obj,tim,tidx,abr,lag) assumes that
%       the responses were pre-aligned with pre-alignment lags, lag, generated 
%       using preproc.m 
%
%       [wabr,tstar] = nlcurvereg(obj,tim,tidx,abr,lag/[],tgt) non-linearly
%       alings the ABRs contained in abr with a predefined structural  
%       target response, tgt, which should contain W elements. If lag = [],
%       the pre-alignment lags are assumed to be zero. 
%
%       [wabr,tstar] = nlcurvereg(obj,tim,tidx,abr,lag/[],[],TGT) if abr is a 
%       cell array with C cells and TGT is a scalar integer between 1 and
%       C, the warping target will be based on the abverage responses for the 
%       TGTth condition (i.e., contained in abr{TGT}). If lag = [],
%       the pre-alignment lags are assumed to be zero. 

    %% Read input.
    vars = {'lag' 'tgt' 'TGT'};
    for I = 1:nargin-4
        eval(sprintf('%s = varargin{%d};',vars{I},I))
    end
    
    %% Create normalized time axis.
    T = numel(tim);
    Dtim = mean(diff(tim));
    
    t = tim(tidx); 
    wrng = [round(t(1)) round(t(end))]
    W = numel(t);
    
    tn = (t-t(1))/(t(end)-t(1)); 

    win = ones(size(t));
    idx = find(and(t>=t(1),t<=t(1)+1));
    win(idx) = sin(pi/2*(t(idx)-t(idx(1)))/(t(idx(end))-t(idx(1))));    
    idx = find(and(t>=t(end)-1,t<=t(end)));
    win(idx) = cos(pi/2*(t(idx)-t(idx(1)))/(t(idx(end))-t(idx(1))));  

    %% Create time warping functions.
    if iscell(abr)
        C = numel(abr);
        N = size(abr{1},1);
        abr = reshape(abr,[1 C]);        
        obj = cell2mat(obj');
    else
        C = 1;
        N = size(abr,1);
    end
    
    if ~exist('lag','var')||isempty(lag)
        if C>1
            lag = zeros(1,C);
        else
            lag = zeros(1,N);
        end
    end
    
    pars = struct('K',diff(wrng),'O',2,'C',5,'LAMBDA',0);
    COST = 'psdd';
    tstar = zeros(N*C,W);
    if ~exist('tgt','var')||isempty(tgt)
        if ~exist('TGT','var')||isempty(TGT)
            NIts = 2;
            tmp = cell(1,NIts);
            for I = 1:NIts, tmp{I} = zeros(N*C,W); end                    
            wobj = obj;
            for I = 1:NIts
                tgt = mean(wobj);
                for II = 1:N*C
                    fprintf(1,'\nit = %d/2 cond = %d/%d subj = %d/%d',I,ceil(II/N),C,II-(ceil(II/N)-1)*N,N);                           
                    tic, [tmp{I}(II,:),wobj(II,:)] = freg(tn,wobj(II,:),tgt,pars,COST); toc
                    if I==2
                        tstar(II,:) = interp1(tn,tmp{I-1}(II,:),tmp{I}(II,:));
                    end   
                end
            end
        else
            NIts = 2;
            tmp = cell(1,NIts);
            for I = 1:NIts, tmp{I} = zeros(N*C,W); end                    
            wobj = obj;
            for I = 1:NIts
                tgt = mean(wobj((TGT-1)*N+1:TGT*N,:));
                for II = 1:N*C
                    fprintf(1,'\nit = %d/2 cond = %d/%d subj = %d/%d',I,ceil(II/N),C,II-(ceil(II/N)-1)*N,N);                           
                    tic, [tmp{I}(II,:),wobj(II,:)] = freg(tn,wobj(II,:),tgt,pars,COST); toc
                    if I==2
                        tstar(II,:) = interp1(tn,tmp{I-1}(II,:),tmp{I}(II,:));
                    end   
                end
            end
        end
    else
        if numel(tgt)==T, tgt = tgt(tidx);  end
        tgt = tgt.*win;
        tgt = tgt-mean(tgt);
        tgt = tgt/sqrt(mean(tgt.^2));
        for I = 1:N*C
            fprintf(1,'cond = %d/%d subj = %d/%d',ceil(I/N),C,I-(ceil(I/N)-1)*N,N);                           
            tic, tstar(I,:) = freg(tn,obj(I,:),tgt,pars,COST); toc
        end
    end    
    
    %% Create warped ABRs.
    if iscell(abr)
        tstar = mat2cell(tstar,ones(1,C)*N); 
        wabr = cellfun(@(x,y) circshift(x,-round(y/Dtim),2),abr,num2cell(lag),'UniformOutput',false);
        for I = 1:C
            for II = 1:N
                wabr{I}(II,tidx) = interp1(tn,abr{I}(II,tidx),tstar{I}(II,:));
            end
        end
    else    
        wabr = abr;
        for I = 1:N
            wabr(I,:) = circshift(wabr(I,:),-round(lag(I)/Dtim));
            wabr(I,tidx) = interp1(tn,abr(I,tidx),tstar(I,:));
        end        
    end
end

%% freg - curve registration.
function [tstar,wobj] = freg(tn,obj,tgt,pars,COST)

    Dtn = mean(diff(tn));
    bO = fbspline(tn,pars.K,pars.O);
    
    options = optimoptions(@fmincon,'MaxFunEvals',2e4,'MaxIter',1e4);
    c = fmincon(@(x)[1 pars.LAMBDA]*fcost(x,tn,obj,tgt,bO,COST),...
        zeros(1,pars.K+pars.O-1),[],[],[],[],-pars.C*ones(1,pars.K+pars.O-1),...
        pars.C*ones(1,pars.K+pars.O-1),[],options);     

    w = sum(bO.*repmat(c',1,length(tn)));
    tstar = cumsum(exp(w)*Dtn); 
    tstar = (tstar-tstar(1))/(tstar(end)-tstar(1));
    
    wobj = interp1(tn,obj,tstar); 
end

%% fcost - cost function.
function Y = fcost(x,tn,obj,tgt,bO,COST)

    Dtn = mean(diff(tn));
    
    w = sum(bO.*repmat(x',1,length(tn)));
    
    tstar = cumsum(exp(w)*Dtn); 
    tstar = (tstar-tstar(1))/(tstar(end)-tstar(1));
    wobj = interp1(tn,obj,tstar);

    Y = zeros(2,1);
    switch(COST)
        case('psd')
        Y(1) = sum((tgt-wobj).^2*Dtn);        
        case('psdd')
        Y(1) = sum((diff(tgt)/Dtn-diff(wobj)/Dtn).^2*Dtn);        
        case('pmc')
        Y(1) = -corr(tgt',wobj');
        case('pmcd')
        Y(1) = -corr(diff(tgt')/Dtn,diff(wobj')/Dtn);            
    end        
    Y(2) = sum((diff(log(diff(tstar)/Dtn))/Dtn).^2*Dtn);
end

%% fbspline - B-spline basis functions.
function bO = fbspline(tn,K,O)

    [~,brk] = min(abs(tn'-(0:1/K:1)));    
    pp = bspline(0:1/O:1);

    bO = zeros(K+O-1,length(tn));
    for I = 1:K+O-1
        for II = 0:O-1
            if and(I-II>0,I-II<=K)
                bO(I,brk(I-II):brk(I-II+1)) = polyval(pp.coefs(O-II,:),(tn(brk(I-II):brk(I-II+1))-tn(brk(I-II)))/(O*(tn(brk(I-II+1))-tn(brk(I-II)))));
            end
        end
    end   
end      

