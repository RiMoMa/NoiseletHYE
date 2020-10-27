clear
load('TestData.mat')
Topic=10;
Epsilon=1e-4;
% fit a plsa model from a given term-doc matrix
tic
[prob_term_topic, prob_topic_doc, lls] = plsa(full(termDocMatrix), Topic, Epsilon);
toc

% plot the log-likelihood per iteration
figure;
plot(lls(2:end));
xlabel('Iteration');
ylabel('log-likehood');
