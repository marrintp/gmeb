function [Y] = diffusionKernel(distMat,sigmaK,alpha,d)
    D = distMat;
    K = exp(-(D/sigmaK).^2);
    p = sum(K);
    p = p(:);
    K1 = K./((p*p').^alpha);
    v = sqrt(sum(K1));
    v = v(:);
    A = K1./(v*v');
    if sigmaK >= 0.5
        thre = 1e-7;  
        M = max(max(A));
        A = sparse(A.*double(A>thre*M));
        [U,~,~] = svds(A,d+1);   %Sparse version.
        U = U./(U(:,1)*ones(1,d+1));
    else
        [U,~,~] = svd(A,0);   %Full version.
        U = U./(U(:,1)*ones(1,size(U,1)));
    end
    Y = U(:,2:d+1);
end