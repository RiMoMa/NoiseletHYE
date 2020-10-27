clear
load('TestData.mat')
Topic=10;
%Epsilon=1e-5;

% fit a plsa model from a given term-doc matrix
tic
[prob_term_topic, prob_topic_doc, lls] = plsaTEM(full(termDocMatrix), Topic, 0.95, 5e-3, 5e-4);
toc

% plot the log-likelihood per iteration
figure;
plot(lls(2:end));
xlabel('Iteration');
ylabel('log-likehood');
