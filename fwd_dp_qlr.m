A=[0.9974 0.0539; -0.1078 1.1591];
B=[0.0013; 0.0539];
H=0;
Q=[0.25 0.00;0.00 0.05];
R=0.05;
N=200;
P=cell(N,1);
P{1}=H;
X=[];
X(:,1)=[2;1];
F=[];
U=[];

for k = 1:N-1

F(N-k,:) = -inv(R+transpose(B)*P{k}*B).*(transpose(B)*P{k}*A);

P{k+1} = transpose(A+B*F(N-k,:))*P{k}*(A+B*F(N-k,:)) + transpose(F(N-k,:))*R*F(N-k,:) + Q;
;

end

for k = 1:N-1

U(k,1) = F(k,:)*X(:,k);
X(:,k+1) = A*X(:,k)+B*U(k,1);

end

U(N,1)=0;
figure
plot(1:size(X,2),X(1,:),'r-',1:size(X,2),X(2,:),'r--',1:size(X,2),U,'b*');
figure
plot(1:size(F,1),F(:,1),'r-',1:size(F,1),F(:,2),'b-');