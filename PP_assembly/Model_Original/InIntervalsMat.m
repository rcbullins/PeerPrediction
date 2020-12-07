function out = InIntervalsMat(ts,X,epochs)
out = cell(size(epochs,1),1);
eps1 = MergeEpochs2(epochs);
kp = (ismember(epochs,eps1,'rows'));
kp1 = ~kp;
kp = find(kp);
[eps2,b] = sortrows(epochs(kp,:));

kp = kp(b);
[status, interval] = InIntervals(ts,eps2);
X1 = X(status);


n1 = histoc(interval(interval>0),1:size(eps2,1));
out(kp) = mat2cell( X1,n1);


out(kp1) = cellfun(@(a) X(ts>=a(1) & ts<a(2)),num2cell(epochs(kp1,:),2),'uni',0);
end