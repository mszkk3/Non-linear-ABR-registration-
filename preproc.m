function [obj,tidx,varargout] = preproc(abr,tim,wrng,varargin)
%preproc ABR preprocessing
%
%   Syntax:
%       [obj,tidx] = preproc(abr,tim,wrng)
%       [obj,tidx,lag] = preproc(abr,tim,wrng,tgt)
%       [obj,tidx,lag] = preproc(abr,tim,wrng,[],TGT)
%       [obj,tidx,lag] = preproc(abr,tim,wrng,tgt/[],TGT/[],lag0)
%       [obj,tidx,lag] = preproc(abr,tim,wrng,tgt/[],TGT/[],[],srng)
%       [obj,tidx,lag] = preproc(abr,tim,wrng,tgt/[],TGT/[],lag0,srng)
%
%   Description:
%       [obj,tidx] = preproc(abr,tim,wrng) preprocesses the ABRs contained in abr by 
%       root-mean-square amplitude normalizing each individual response.
%       abr ==> is either an array with dimensions N (number of subjects) by T
%       (number of time points), or a cell array with C cells, where C is
%       the number of conditions, and each cell has dimensions N by T.
%       tim ==> is an array of the T time samples in milliseconds. 
%       wrng is a two-element vector specifying the warping time range
%       (e.g, [0 12]).
%       obj ==> is an array/cell-array of normalized responses, ready for
%       warping.
%       tidx ==> indeex of time samples contained within the warping time
%       range, wrng.
%       
%       [obj,tidx,lag] = preproc(abr,tim,wrng,tgt) will additionally pre-align the
%       responses in abr to a predefined target response, tgt. 
%       tgt ==> must contain either T or W elements, where W is the length of 
%       the warping time range in samples (tidx). 
%       lag ==> is the N-or C-element vector of pre-alignement lags (if abr is a cell 
%       array, the pre-alignment lags will be based on the average responses across 
%       subjects for each condition).
%
%       [obj,tidx,lag] = preproc(abr,tim,wrng,[],TGT) if abr is a cell array with C cells  
%       and TGT is a scalar integer between 1 and C, the pre-aligment will be performed 
%       using the average response for the TGTth condition (i.e., contained in abr{TGT})
%       as pre-alignment target. 

%       [obj,tidx,lag] = preproc(abr,tim,wrng,tgt/[],TGT/[],lag0) limits the search for
%       pre-alignment lags to the +/-2-ms time range around lag0.
%       lag0 ==> is a N-or C-element vector of initial lags. By default
%       (if lag0 is either empty or not defined), lag0 is set to zero. 
%
%       [obj,tidx,lag] = preproc(abr,tim,wrng,tgt/[],TGT/[],[],srng) limits the search for 
%       pre-alignment lags to the specified time range around zero. 
%       srng ==> is a two-element vector specifying the lower and upper bounds
%       of the range in milliseconds (e.g., [0 2]). By default (if srng is
%       not defined), srng is set to [-2 2].
%
%       [obj,tidx,lag] = preproc(abr,tim,wrng,tgt/[],TGT/[],lag0,srng) limits the search for 
%       pre-alignment lags to the specified time range around lag0. 

    %% Read input.
    vars = {'tgt' 'TGT' 'lag0' 'srng'};
    for I = 1:nargin-3
        eval(sprintf('%s = varargin{%d};',vars{I},I))
    end

    %% Create warping time window.
    T = numel(tim);
    Dtim = mean(diff(tim));
    
    wrng = wrng(1)+[0 round(diff(wrng))];
    tidx = find(and(tim>=wrng(1),tim<=wrng(2)));     

    t = tim(tidx);
    W = numel(t);
    win = ones(size(t));
    idx = find(and(t>=t(1),t<=t(1)+1));
    win(idx) = sin(pi/2*(t(idx)-t(idx(1)))/(t(idx(end))-t(idx(1))));    
    idx = find(and(t>=t(end)-1,t<=t(end)));
    win(idx) = cos(pi/2*(t(idx)-t(idx(1)))/(t(idx(end))-t(idx(1))));  
        
    %% Find pre-alignment lags.
    if iscell(abr)
        C = numel(abr);
        N = size(abr{1},1);
        abr = reshape(abr,[1 C]);
        lag = zeros(1,C); 
    else
        C = 1;
        N = size(abr,1);
        lag = zeros(1,N); 
    end
    
    if exist('tgt','var')&&~isempty(tgt)
        if numel(tgt)==T, tgt = tgt(tidx);  end
        tgt = reshape(tgt,[1 W]);
    else
        tgt = [];
    end
    
    if exist('TGT','var')&&~isempty(TGT)
        tgt = mean(abr{TGT}(:,tidx));
    end

    if ~isempty(tgt)
        if iscell(abr)
            tmp = cell2mat(cellfun(@(x) mean(x),abr','UniformOutput',false));
            if ~exist('lag0','var')||isempty(lag0), lag0 = zeros(1,C); end                
        else
            tmp = abr;
            if ~exist('lag0','var')||isempty(lag0), lag0 = zeros(1,N); end
        end
        if ~exist('srng','var'), srng = [-2 2]; end

        for I = 1:size(tmp,1)
            lags = lag0(I)+srng(1):Dtim:lag0(I)+srng(2);
            y = cell2mat(arrayfun(@(x) circshift(tmp(I,:),-round(x/Dtim))',lags,'UniformOutput',false));
            
            [~,IDX] = max(corr((tgt.*win)',y(tidx,:)));
            lag(I) = lags(IDX);
        end
        varargout = {lag};
    end
    
    %% Preprocess responses.        
    if iscell(abr)
        obj = cellfun(@(x,y) circshift(x,-round(y/Dtim),2),abr,num2cell(lag),'UniformOutput',false);
        obj = cellfun(@(x) x(:,tidx).*repmat(win,N,1),obj,'UniformOutput',false);
        obj = cellfun(@(x) x-mean(x,2),obj,'UniformOutput',false);
        weight = cellfun(@(x) 1./sqrt(mean(x.^2,2)),obj,'UniformOutput',false); 
        obj = cellfun(@(x,y) x.*repmat(y,1,W),obj,weight,'UniformOutput',false); 
    else
        obj = cell2mat(cellfun(@(x,y) circshift(x,-round(y/Dtim),2),num2cell(abr,2),num2cell(lag)','UniformOutput',false));
        obj = obj(:,tidx);
        obj = obj-mean(obj,2);
        weight = 1./sqrt(mean(obj.^2,2));
        obj = obj.*repmat(weight,1,W);
    end
end

