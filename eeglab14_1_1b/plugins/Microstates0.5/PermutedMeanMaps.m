function [MeanMap,SortedMaps,OldMapFit] = PermutedMeanMaps(in,RespectPolarity, Montage)

    if nargin < 3
        debug = false;
    else
        debug = true;
    end
    
    [nSubjects,nMaps,~] = size(in);

    progStrArray = '/-\|';
    fprintf(1,'Permuting %i maps of %i subjects ',nMaps,nSubjects);

    in = NormDim(in,3);
    if nargin < 2;  RespectPolarity = false;    end

    MeanMap = in(1,:,:);
    SortedMaps = in;
    OldSortedMaps = SortedMaps;
    OldMapFit = -inf;
    WorkToBeDone = true;
    cnt = 0;
    
    NewOrder = repmat(1:nMaps,nSubjects,1);
    OldOrder = NewOrder;
    OldOldOrder = OldOrder;
            
    while WorkToBeDone
        cnt = cnt + 1;
        fprintf(1,'\b%s',progStrArray(mod(cnt-1, 4)+1));
        % See how the prototype fits
        MapFit = GetMapSeriesOverallFit(SortedMaps,MeanMap,RespectPolarity);         
        
        if debug == true
            figure(1000);
            spidx = 1;
            X = cell2mat({Montage.X});
            Y = cell2mat({Montage.Y});
            Z = cell2mat({Montage.Z});
            for s = 1:nSubjects
                for c = 1:nMaps
                    subplot(nSubjects+1,nMaps,spidx);
                    dspQMap(squeeze(SortedMaps(s,c,:)),[X;Y;Z],'Resolution',5);
                    spidx = spidx+1;
                end
            end
            for c = 1:nMaps
                subplot(nSubjects+1,nMaps,spidx);
                dspQMap(squeeze(MeanMap(1,c,:)),[X;Y;Z],'Resolution',5);
                spidx = spidx+1;
            end
        end
        
        MeanMapFit = mean(MapFit);
            
        % Find the order of misfit
        [~,Idx] = sort(MapFit(:),'ascend');
        WorkToBeDone = false;
        for i = 1:numel(Idx)
            if (nMaps < 7) || (license('test','optimization_toolbox') == false) % Full permutations for small n or absent optimzation toolbox
                [SwappedMaps,Order] = SwapMaps(SortedMaps(Idx(i),:,:),MeanMap,RespectPolarity);
            else        % linear prgramming for larger problems
                [SwappedMaps,Order] = SwapMaps2(SortedMaps(Idx(i),:,:),MeanMap,RespectPolarity);
            end
            if ~isempty(SwappedMaps)
                if debug == true
                    Idx(i)
                end
                WorkToBeDone = true;
                SortedMaps(Idx(i),:,:) = SwappedMaps;
                NewOrder(Idx(i),:) = Order;
                for k = 1:nMaps
                     data = squeeze(SortedMaps(:,k,:));
                     [pc1,~] = eigs(data'*data,1);
                     MeanMap(1,k,:) = NormDimL2(pc1',2);
                end
                break;
            else
                if debug == true
                    fprintf(1,'No change in case %i\n',Idx(i));
                end
            end
            
        end
       
        % Catch the rate cases where the linintprog swaps back and forth
        if all(OldOldOrder(:) == NewOrder(:))
            disp('flip back found');
            SortedMaps = OldSortedMaps;
            for k = 1:nMaps
                 data = squeeze(SortedMaps(:,k,:));
                 [pc1,~] = eigs(data'*data,1);
                 MeanMap(1,k,:) = NormDimL2(pc1',2);
            end
            break;
        end

        OldOldOrder = OldOrder;
        OldOrder = NewOrder;
        
        OldSortedMaps = SortedMaps;
        OldMapFit     = MeanMapFit;
    end
    
    if debug == true
        disp('Done');
        pause
    end
    for s = 1:size(SortedMaps,1)
        SubMaps = squeeze(SortedMaps(s,:,:));
        for k = 1:size(SortedMaps,2)
            if  SubMaps(k,:) *  squeeze(MeanMap(1,k,:)) < 0
                SortedMaps(s,k,:) = -SortedMaps(s,k,:);
            end
        end
    end
    
    fprintf(1,'\b, mean fit: %f\n',mean(MapFit(:)));
 
    MeanMap = squeeze(MeanMap); 
    covm = MeanMap'*MeanMap;
    [v,d] = eigs(covm,1);
    sgn = sign(MeanMap * v);
    MeanMap = MeanMap .* repmat(sgn,1,size(MeanMap,2));
        
end


