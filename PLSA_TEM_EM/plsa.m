function [prob_term_topic, prob_topic_doc,lls] = plsa(termDocMatrix, numTopic, epsilon)

% Fit a plsa model from a given term-document matrix with ordinal EM
% algortihm
% This code was inspired by the following code (https://github.com/lizhangzhan/plsa)
% AUTHOR: Ikko Kimura, Osaka University, 2016/04/02
% USAGE: [prob_term_topic, prob_topic_doc,lls] = plsa(termDocMatrix, numTopic, epsilon)
% epsilon: criteria for stoppping EM algorithm (default: 1e-5)
% As for pLSA please refer to the following paper
% Thomas Hofmann, 2001, Machine Learning, 42, 177-196, Unsupervised Learning by
% Probabilistic Latent Semantic Analysis

if nargin < 3
    epsilon=1e-5;
end

[numTerm, numDoc] = size(termDocMatrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!! INITIALIZATIONS !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
prob_term_topic = rand(numTerm, numTopic); % p(term | topic)
prob_term_topic = bsxfun(@rdivide,prob_term_topic,sum(prob_term_topic));

prob_topic_doc = rand(numTopic, numDoc);   % p(doc | topic)
prob_topic_doc = bsxfun(@rdivide,prob_topic_doc,sum(prob_topic_doc));

prob_topic_term_doc=zeros(numTerm, numDoc,numTopic);

prob_term_doc=prob_term_topic * prob_topic_doc;

lls =sum(sum(termDocMatrix.* log(prob_term_doc+eps)));  % maximum log-likelihood estimations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!! ITERATION STARTS !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
done=0;

while ~done; i = i+ 1; 
	%%% E-step...
    [prob_topic_term_doc]=Estep(prob_topic_doc,prob_term_topic,prob_term_doc,prob_topic_term_doc);
    %%% M-step...
    [prob_topic_doc,prob_term_topic]=Mstep(termDocMatrix,prob_topic_term_doc);
    
    %%% calculate likelihood and update p(term, doc) 
    prob_term_doc  = prob_term_topic*prob_topic_doc;
    ll=sum(sum(termDocMatrix.* log(prob_term_doc+eps)));
	
    display(sprintf('Iteration %d Loglikelihood %d',i,ll));
    lls= [lls;ll];
    
    rel_ch = (lls(i) - lls(i-1))/abs(lls(i-1)); 
    if rel_ch< epsilon && i> 5 
        done=1;
        display(sprintf('Convergence was made after %d steps..',i))
    end
end

end

function [prob_topic_term_doc]=Estep(prob_topic_doc,prob_term_topic,prob_term_doc,prob_topic_term_doc)
%perform E step
for z = 1:size(prob_topic_term_doc,3) %% Hmm...
prob_topic_term_doc(:, :, z) = bsxfun(@times,prob_topic_doc(z,:),prob_term_topic(:,z))./(prob_term_doc+eps);
end
%prob_topic_term_doc(isnan(prob_topic_term_doc))=0;
end

function [prob_topic_doc,prob_term_topic]=Mstep(termDocMatrix,prob_topic_term_doc)
% perform M step
g=bsxfun(@times,termDocMatrix,prob_topic_term_doc);
prob_topic_doc=permute(sum(g,1),[3 2 1]);
prob_term_topic=permute(sum(g,2),[1 3 2]);
% Normalize the data
prob_topic_doc = bsxfun(@rdivide,prob_topic_doc,sum(prob_topic_doc));
prob_term_topic = bsxfun(@rdivide,prob_term_topic,sum(prob_term_topic));
end