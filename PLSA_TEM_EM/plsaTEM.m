function [prob_term_topic, prob_topic_doc,lls] = plsaTEM(termDocMatrix, numTopic, nu, epsilon1, epsilon2)

% Fit a plsa model from a given term-document matrix with Tempered EM
% algortihm
% This code was inspired by the following code (https://github.com/lizhangzhan/plsa)
% AUTHOR: Ikko Kimura, Osaka University, 2016/04/02
% 2016/04/09: change the stopping critria to epsion2/t in every bÅ©b*nu steps
% 2016/04/10: made the code much faster

% USAGE: [prob_term_topic, prob_topic_doc,lls] = plsaTEM(termDocMatrix, numTopic, nu, epsilon1, epsilon2)
% nu: parameter for betaÅ®nu*beta (default: 0.95)
% epsilon1: criteria for early stopping EM (default: 1e-3)
% epsilon2: criteria for stoppping EM algorithm (default: 5e-3)
% As for Tempered EM and pLSA please refer to the following paper
% Thomas Hofmann, 2001, Machine Learning, 42, 177-196, Unsupervised Learning by
% Probabilistic Latent Semantic Analysis

if nargin < 5
    epsilon2=5e-3;
    if nargin <4
        epsilon1=1e-3;
    if nargin <3
        nu=0.95;
    end
    end
end

[numTerm, numDoc] = size(termDocMatrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!! INITIALIZATIONS !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
prob_term_topic = rand(numTerm, numTopic); % p(term | topic)
prob_term_topic = prob_term_topic./repmat(sum(prob_term_topic),[numTerm,1]);

prob_topic_doc = rand(numTopic, numDoc);   % p(doc | topic)
prob_topic_doc = prob_topic_doc./repmat(sum(prob_topic_doc),[numTopic,1]);

prob_topic_term_doc=zeros(numTerm, numDoc,numTopic);

prob_term_doc=prob_term_topic * prob_topic_doc;

lls =sum(sum(termDocMatrix.* log(prob_term_doc+eps)));  % maximum log-likelihood estimations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!! ITERATION STARTS !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
done=0;

%%% FIRST PERFORM EM WITH EARLY STOPPING...
display('PERFORM EM WITH EARLY STOPPING...')
while ~done; i = i+ 1; 
	%%% E-step...
	[prob_topic_term_doc]=Estep1(prob_topic_doc,prob_term_topic,prob_term_doc,prob_topic_term_doc);
	%%% M-step...
	[prob_topic_doc,prob_term_topic]=Mstep(termDocMatrix,prob_topic_term_doc);
    
    %%% calculate likelihood and update p(term, doc) 
    prob_term_doc  = prob_term_topic*prob_topic_doc;
    ll=sum(sum(termDocMatrix.* log(prob_term_doc+eps)));
	
    display(sprintf('Iteration1 %d Loglikelihood %d',i,ll));
    lls= [lls;ll];
    
    rel_ch = (lls(i) - lls(i-1))/abs(lls(i-1)); 
    if rel_ch< epsilon1
        done=1;
        display(sprintf('Convergence was made after %d steps..',i))
    end
end

%%% DO THE SECOND STEPS...
Done=0;
b=1; % beta
t=0;
while ~Done; 
    b=b*nu;
    t=t+1;
    display(sprintf('Now the beta is %d',b));
    k=0; done=0;
while ~done; 
    i = i+ 1; k = k+1;
    
	%%% E-step...
	[prob_topic_term_doc]=Estep2(prob_topic_doc,prob_term_topic,prob_topic_term_doc,b);
	%%% M-step...
	[prob_topic_doc,prob_term_topic]=Mstep(termDocMatrix,prob_topic_term_doc);
    
    %%% calculate likelihood and update p(term, doc) 
    prob_term_doc  = prob_term_topic*prob_topic_doc;
    ll=sum(sum(termDocMatrix.* log(prob_term_doc+eps)));
    display(sprintf('Iteration%d %d Loglikelihood %d',t+1,k,ll));
    
    lls= [lls;ll];
    
    rel_ch = (lls(i) - lls(i-1))/abs(lls(i-1)); 
    if rel_ch< epsilon2/t
        if k==1
            Done=1;
        end
        done=1;
    end
end
end


end

function [prob_topic_term_doc]=Estep1(prob_topic_doc,prob_term_topic,prob_term_doc,prob_topic_term_doc)
%perform ordinal E step
for z = 1:size(prob_topic_term_doc,3) %% Hmm...
prob_topic_term_doc(:, :, z) = bsxfun(@times,prob_topic_doc(z,:),prob_term_topic(:,z))./ (prob_term_doc+eps);
end
%prob_topic_term_doc(isnan(prob_topic_term_doc))=0;
end

function [prob_topic_term_doc]=Estep2(prob_topic_doc,prob_term_topic,prob_topic_term_doc,b)
%perform tempered E step
for z = 1:size(prob_topic_term_doc,3) %% Hmm...
prob_topic_term_doc(:, :, z) =(bsxfun(@times,prob_topic_doc(z,:),prob_term_topic(:,z))).^b;
end
%prob_topic_term_doc(isnan(prob_topic_term_doc))=0;
prob_topic_term_doc = bsxfun(@rdivide,prob_topic_term_doc,sum(prob_topic_term_doc,3)+eps);
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