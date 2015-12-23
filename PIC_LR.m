function [test_rmses,timelist,iterlist,objvallist]=PIC_LR(m,n,GXobs,GXtest, maxiter, maxk,b)

gamma=10000;lambda=0.1;
R = GXobs;
kk=(maxk/b)*ones(1,b);
Udm={};Vdm={};
for i=1:b
    Udm{i}=zeros(m,kk(i));
    Vdm=[Vdm;zeros(kk(i),n)];
end
iterlist=[0];
countiter=0;
test_rmses=[1];
eto_test=1e-4;
etol=1e-2;
loop_num=800;
for loop=1:loop_num
    etol=max(etol*0.9,1e-4);
    test_rmses_old=10000;
    %lambda1=max(eta*lambda,1e-2);gamma1=gamma;
    if loop==1200
        maxiter=100*b;
    end
    %[Udm1,Vdm1,R1,test_rmses1,timelist1,iterlist1,countiter1,objvallist1]=cdmnnMC(m,n,GXobs,Udm,Vdm,R,lambda1,gamma1,GXtest,test_rmses,iterlist,countiter, maxiter, maxk,b,etol);
    lambda2=lambda;gamma2=max(0.95*gamma,1.05);
    [Udm2,Vdm2,R2,test_rmses2,timelist2,iterlist2,countiter2,objvallist2]=cdmnnMC(m,n,GXobs,Udm,Vdm,R,lambda2,gamma2,GXtest,test_rmses,iterlist,countiter, maxiter, maxk,b,etol);
    if 0%objvallist1(end)<objvallist2(end)
        Udm=Udm1;Vdm=Vdm1;R=R1;lambda=lambda1;gamma=gamma1;
        test_rmses=test_rmses1;timelist=timelist1;iterlist=iterlist1;countiter=countiter1;
        objvallist=objvallist1;
    else
        Udm=Udm2;Vdm=Vdm2;R=R2;lambda=lambda2;gamma=gamma2;
        test_rmses=test_rmses2;timelist=timelist2;iterlist=iterlist2;countiter=countiter2;
        objvallist=objvallist2;
    end
    if abs(test_rmses(end)-test_rmses_old)/max(test_rmses_old,1)<eto_test
        return;
    end
    test_rmses_old=test_rmses(end);
end
