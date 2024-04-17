% Removes all rows and columns that only contain nans
function A=nanremove(A)

A(:,all(isnan(A),1))=[];
A(all(isnan(A),2),:)=[];

