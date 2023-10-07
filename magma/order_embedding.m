// hash function walk from given order, with steps of degree L defined by vector steps (this vector must contain elements from 0 to L-1)
HashWalk:=function(order,L,steps)
	// quaternion algebra
	B<i,j,k>:=Parent((Basis(order))[1]);
	p:=-Integers()!Rationals()!(i^2);
	q:=-Integers()!Rationals()!(j^2);

	// preparing for a random walk
	myprecision:=Ceiling(10*Log(L,p));
	BL,B2BL,Q2QL:=pMatrixRing(order,L: Precision:=myprecision);
	//BL,B2BL,Q2QL:=pMatrixRing(order,L);
	BL2B:=Inverse(B2BL);
	BL:=Codomain(B2BL);
	QL:=Codomain(Q2QL); 
	ZL:=IntegerRing(QL);
	vecquat:=[];
	for ind in [1..L] do
		vecquat:=vecquat cat [BL2B(Matrix(ZL,2,2,[L,0,ind-1,1]))];
	end for;
	vecquat:=vecquat cat [BL2B(Matrix(ZL,2,2,[1,0,0,L]))];
	//"norm of vecquat elements have L valuation equal to",[Valuation(Norm(quat),L): quat in vecquat];

	// hash function walk
	length:=#steps;
	alpha:=1;N:=1;
	for ind in [1..length] do
		step:=steps[ind];
		N:=N*L;
		alpha:=alpha*vecquat[step+1];
		//alpha:=(Integers()!alpha[1] mod N) + (Integers()!alpha[2] mod N)*i + (Integers()!alpha[3] mod N)*j + (Integers()!alpha[4] mod N)*k;
	end for;
	I:=LeftIdeal(order,[N,alpha]);
	order:=RightOrder(I);

	return order;
end function;



// random walk from given order, with N random steps of degree L
RandomWalk:=function(order,L,length)
	order:=HashWalk(order,L,[Random(L-1): ind in [1..length]]);
	return order;
end function;

RunMeOnce := function(p, steps)
	// p is prime, steps was bit length of p
	Q:=Rationals();
	B<i,j,k>:=QuaternionAlgebra<Q|-1,-p>;
	O:=QuaternionOrder([1,i,j,k]);
	O1:=MaximalOrder(O);
	O2:=RandomWalk(O1,2,steps);
	s1:=Basis(O2)[2];
	s2:=Basis(O2)[3];
	s3:=Basis(O2)[4];
	Trace(s1);
	Trace(s2);
	Trace(s3);
	S:=[];
	S[1]:=s1;
	S[2]:=s2;
	S[3]:=s3;
	M:=Matrix(Q,3,3, [Trace(Conjugate(S[i])*S[j]):i,j in [1..3]]);
	N:=LLLGram(M);
	return 1/2*QuadraticForm(N);
end function;

RunMeNTimes := function(p, steps, n)
	res := [];
	for i in [1..n] do
		form := RunMeOnce(p, steps);
		res[i] := form;
	end for;
	return res;
end function;
