function [P] = computePropensities(k,X)
	P=zeros(size(X,1),6);
	P(:,1) = k(:,2).*X(:,1);
	P(:,2) = k(:,1).*X(:,2);
	P(:,3) = k(:,3).*X(:,1);
	P(:,4) = k(:,4).*X(:,3);
	P(:,5) = k(:,5).*X(:,3);
	P(:,6) = k(:,6).*X(:,4);

end