function [SortedMaps,SortOrder, Communality, Polarity] = ArrangeMapsBasedOnMean(in, MeanMap,RespectPolarity)

    [nSubjects,nMaps,nChannels] = size(in);

    Communality = nan(nSubjects,nMaps);
    Polarity    = ones(nSubjects,nMaps);
    fprintf(1,'\nSorting %i maps of %i subjects. ',nMaps,nSubjects);

    in = NormDim(in,3);
  
    if nargin < 3
        RespectPolarity = false;
    end

    SortedMaps = in;
    ExtMeanMap(1,:,:) = MeanMap;
    SortOrder = repmat(1:nMaps,nSubjects,1);
    
    for n = 1:nSubjects
		MapsToSort = in(n,:,:);

        if (nMaps < 7) || (license('test','optimization_toolbox') == false) % Full permutations for small n or absent optimzation toolbox
            [SwappedMaps,Assignment, pol] = SwapMaps(MapsToSort,ExtMeanMap,RespectPolarity);
        else        % linear prgramming for larger problems
            [SwappedMaps,Assignment,pol] = SwapMaps2(MapsToSort,ExtMeanMap,RespectPolarity);
        end
        Polarity(n,:) = pol;
        if ~isempty(SwappedMaps)
            SortedMaps(n,:,:) = SwappedMaps;
            SortOrder(n,:) = Assignment;
        end

        Communality(n,:) = diag(corr(squeeze(SortedMaps(n,:,:))',squeeze(ExtMeanMap)'))';
    end
    if ~RespectPolarity
        Communality = abs(Communality);
    end
end
