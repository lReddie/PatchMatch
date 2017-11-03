function [res, ann, bnn] = search_vote_funcnaive(img, res, niters)%, j, k, resultFolder, sceneName)
% Please replace the nnmex and votemex after you have implemented your own
% searching and voting function for PatchMatch

%% Parameters for patchmatch

%% Go th n iterations
for iter = 1:niters
    %% Searching for NNF
        ann = naive(res, img);
        bnn = naive(img, res);

	%% Voting     
    res = votemex(img, ann, bnn);    

end
end