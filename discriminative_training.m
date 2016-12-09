function classLoss=discriminative_training(X, Y)
    rng(1);
    SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','linear', 'KernelScale',1, 'CVPartition', cvpartition(Y,'kfold', 5));
    %CVSVMModel = crossval(SVMModel);
    classLoss = kfoldLoss(SVMModel,'lossfun',@loss);
end

function l=loss(ground_truth,predictions,W,cost)
assert(length(predictions)==length(ground_truth),'lengths do not match\n');
ground_truth=ground_truth(:,1);
[~,idx]=max(predictions,[],2);
predictions=(idx==1);
TP=sum((predictions==ground_truth).*(predictions==1));
TN=sum((predictions==ground_truth).*(predictions==0));
FP=sum((predictions~=ground_truth).*(predictions==1));
FN=sum((predictions~=ground_truth).*(predictions==0));

precision=TP/(TP+FP);recall=TP/(TP+FN);
if precision==0 || recall==0
    F1=0;
else
    F1=2*(precision*recall)/(precision+recall);
end
%l=F1;
%l=TP/(TP+FN); % recall
%l=TP/(TP+FP);  % precision
l=(TP+TN)/(FP+FN+TP+TN); % accuracy
end
