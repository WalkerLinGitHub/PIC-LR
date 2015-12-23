function [Udm,Vdm,R,test_rmses,timelist,iterlist,countiter,objvallist]=cdmnnMC(m,n,GXobs,Udm,Vdm,R,lambda,gamma,GXtest,test_rmses,iterlist,countiter, maxiter, maxk,b,etol)

%obj_vals=zeros(maxiter,1);
%GXobs=sparse(GXobs);
%[i_row, j_col, data]=find(GXobs); number_observed=length(data);
[i_row_test, j_col_test, data_test]=find(GXtest); number_observed_test=length(data_test);
timelist=[0];objvallist=[0];
% initial from zero
kk=(maxk/b)*ones(1,b);

totaltime = 0;
oldV = [];

compRtime=0;
i=0; 
choseni=0;
Tmax=3; Xp=zeros(m,n);
while (i<maxiter) 
	timebegin = cputime;
	i=i+1;

    choseni=mod(choseni,b)+1;
    if choseni==1
        countiter=countiter+1;
    end

    R1 = CompResidual(R, -Udm{choseni}', Vdm{choseni});
    Li=0.65;
    %[u,s,v]=svd(1.0/(Li)*R+Udm{choseni}*Vdm{choseni});%1.1< 1.2 <1.4
    %u=u(:,1:kk{choseni});s=s(1:kk{choseni},1:kk{choseni});v=v(:,1:kk{choseni});
    %fprintf('CDMNN,choseni:%g, R: %g\n',choseni,norm(R,'fro')); 
    [u,s,v] = randomsvd(R/Li, Udm{choseni}, Vdm{choseni}', m, n, kk(choseni), oldV, Tmax);
    %[u,s,v] = svds(R/Li+Udm{choseni}*Vdm{choseni}, kk(choseni));
    
	sing_vals=diag(s); clear s;

    tmp=sing_vals;
    tmp(sing_vals<=lambda)=0;
    ind=find((lambda<sing_vals)&&(sing_vals<=lambda*gamma));
    tmp(ind)=gamma/(gamma-1)*(sing_vals(ind)-lambda/Li);
	%tmp=max(sing_vals-lambda/Li,0);
      
	soft_singvals=tmp(tmp>0);
    
    S=diag(soft_singvals);
	U = u(:,tmp>0); clear u;
	V = v(:,tmp>0); clear v;
    if size(U,2)>kk(choseni)
        U=U(:,1:kk(choseni));S=S(1:kk(choseni),1:kk(choseni));V=V(:,1:kk(choseni));
    end
    oldV = V;
    
    V=S*V';
    
    Udm{choseni}=U;
    Vdm{choseni}=V;
    
    t2=cputime;
    R = CompResidual(R1, Udm{choseni}', Vdm{choseni});
    compRtime=compRtime+cputime-t2;
    totaltime = totaltime + cputime - timebegin;
	
    %fprintf('Iter %g time %g tiker %g choseni %g test_rmse %g\n', i, totaltime, tiker,choseni, test_err);
    
    if choseni==b  
        U0=cell2mat(Udm);
        V0=cell2mat(Vdm);
        xtemp=U0*V0;
        Xp1=xtemp;
        iterlist=[iterlist,countiter];
        
        %Rtemp = CompResidual(GXobs, U0', V0);
        %[~,ssss,~] = svds(xtemp,maxk); train_err = norm(Rtemp,'fro')^2;
        %objval_new = train_err/2 + lambda*sum(diag(ssss));  %lambda
        %objvallist=[objvallist,objval_new];
        
        timelist=[timelist,totaltime];
        
        %%xtemp =  xtemp-GXtest; 
        Rxtemp = CompResidual(GXtest, U0', V0);
        %test_err = sqrt(tttemp_test'*tttemp_test/numel(data_test));
        test_err =  norm(Rxtemp,'fro')/norm(data_test,'fro');
        test_rmses=[test_rmses,test_err];
        %fprintf('Iter %g compRtime %g totaltime %g, train_err %g\n', countiter, compRtime,totaltime,train_err);
        fprintf('Iter %g inner iter %g compRtime %g blocktime %g\n', countiter, i/b,compRtime,totaltime);
        
        if norm(Xp1-Xp,'fro')/norm(Xp,'fro') < etol%test_err < etol 
            %fprintf('Iter %g time %g test_rmse %g\n', i, totaltime, test_err);
            return;
        end
        Xp=Xp1;  
    end
end

