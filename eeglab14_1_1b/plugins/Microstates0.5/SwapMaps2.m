function [Result,Assignment, pol] = SwapMaps2(MapsToSwap,MeanMap,RespectPolarity)
    nMaps = size(MeanMap,2);
    CorMat = corr(squeeze(MapsToSwap)',squeeze(MeanMap)');
    if ~RespectPolarity == true
        CorMat = abs(CorMat);
    end
    DistMat = 1-CorMat(:);
    intcon = 1:numel(DistMat);
    
    % Now we make a weighting matrix that allows us to constrain the thing
    % to one to one assignments
    A = zeros(2*nMaps,numel(DistMat));
    for m = 1:nMaps
        Am1 = zeros(nMaps);
        Am2 = zeros(nMaps);
        Am1(m,:) = 1;
        Am2(:,m) = 1;
        A(m      ,:) = Am1(:);
        A(m+nMaps,:) = Am2(:);
    end
    lb = zeros(numel(DistMat),1);
    ub = ones(numel(DistMat),1);
    b = ones(2*nMaps,1);
    x = intlinprog(DistMat,intcon,[],[],A,b,lb,ub,optimoptions('intlinprog','Display','off'));
    AssignMatrix = reshape(x,nMaps,nMaps);
    [~,Assignment] = max(AssignMatrix == 1);
    
    if all(Assignment == 1:nMaps)
        Result = [];
    else
        Result = MapsToSwap(:,Assignment,:);
        pol = sign(diag(corr(squeeze(Result)',squeeze(MeanMap)')));
    end
end
