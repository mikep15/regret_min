function X=Solve_BQME(A,B,C,D) 
% This Function solves a bilateral matrix quadratic equation 
% of the form AX+XB+XCX+D = 0 for X
% Inputs : Matrices A,B,C, D of appropriate dimensions
%
% Output : The Matrix X - if a solution exists
%
%
% Example
% A =
%      1     2
%      3     5
%      
% B =
%      2     2     1
%      2     3     2
%      2     1     1
% 
% C =
%      1     4
%      1     3
%      2     3
% 
% D =
%      1     5     3
%      2     6     2
% X=Solve_BQME(A,B,C,D)
% 
% X =
%     1.0388   -1.8993    0.1548
%    -0.3030   -0.7949   -0.3434


if ~isreal(A) || ~isreal(B) || ~isreal(C) || ~isreal(D)
    error('Real Valued Matrices Expected');
end

[p1,p2]=size(A);
[p3,p4]=size(B);

if (p1~=p2) || (p3~=p4)
    error('Matrices A and B should be square');
end

if ~isequal(size(C),[p3 p2]) || ~isequal(size(D),[p1 p4])
    error('Matrix Dimensions do not match');
end

X=zeros(p2,p3);
iters=0;
Omega=EvalBiQuad(A,B,C,D,X);
while (iters<100*sum(size(X))) && (norm(Omega)>1e-9) 
    Delta = sylv(A+X*C,B+C*X,-Omega);
    if norm(Delta) > 1e10
        error('Cannot Find Real Valued Solution')
    end
    iters=iters+1;
    X=X+Delta;
    Omega=EvalBiQuad(A,B,C,D,X);
end

end

function Omega=EvalBiQuad(A,B,C,D,X)
Omega=A*X+X*B+X*C*X+D;
end


