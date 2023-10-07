Q:=Rationals();
Generate_Order_Large_Prime_Conductor := function(B,O0,D)
    R<x1,x2,x3,x4> := PolynomialRing(Rationals(), 4); 
    P<X4> := PolynomialRing(GF(D),1);

    M:=Matrix(Q,4,4,[Trace(Conjugate(Basis(O0)[i])*Basis(O0)[j]):i,j in [1..4]]);
    N:=LLLGram(M);
    norm_form:=QuadraticForm(N);

    //guess x1, x2, x3
    repeat 
        x1 := Random(D);
        x2 := Random(D);
        x3 := Random(D);
        f := hom<R -> P | x1,x2,x3,X4>(norm_form);
        if IsIrreducible(f) eq false then
            X4 := GF(D)!-Coefficients(Factorisation(f)[1][1])[2];
        end if;
    until IsIrreducible(f) eq false;
    x4 := Integers()!X4;
    alpha := x1*Basis(O0)[1] + x2*Basis(O0)[2] + x3*Basis(O0)[3] + x4*Basis(O0)[4];
    O2 := RightOrder(LeftIdeal(O0,[alpha,D]));
    I:=LeftIdeal(O0,[alpha,D]);
    return I;
end function;


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
    end for;
    I:=LeftIdeal(order,[N,alpha]);
    order:=RightOrder(I);
    return order,I, alpha;
end function;

RandomWalk:=function(order,L,length)
    order,I, alpha:=HashWalk(order,L,[Random(L-1): ind in [1..length]]);
    return order,I, alpha;
end function;

// Function works for starting order O0, otherwise might have to adapt the multiplication by 2 scaling
Random_NormEqwSol:=function(B, O1, D)
    I :=Generate_Order_Large_Prime_Conductor(B,O1, D);
    // Compute quadratic form
    R<x1,x2,x3,x4> := PolynomialRing(Rationals(), 4);
    M := ChangeRing(ReducedGramMatrix(I)/D, R);

    quadratic_form := (Matrix(R, 1, 4, [x1, x2, x3, x4]) * M * Matrix(R, 4, 1, [x1, x2, x3, x4]))[1][1];
    // Warning!: Compared to ReducedBasis, ReducedGramMatrix seems to multiply all coefficients by 2 (to get integer coefficients?), therefore we divide by 2 to get something consistent with ReducedBasis.
    quadratic_form := quadratic_form / 2;

    // Compute solution for check (note this uses implicitly that I is of norm L^steps already)
    redbasis := ReducedBasis(I);
    coeffs := [];
    for i in [1..4] do
        coeffs := coeffs cat Eltseq(redbasis[i]);
    end for;
    CoeffMat := Matrix(Rationals(),4,4,coeffs);
    sol := Solution(CoeffMat,Matrix(Rationals(),1,4,[D,0,0,0]));

    // Check that solution has norm L^steps
    if Evaluate(quadratic_form, [Integers()!sol[1][j]: j in [1..4]]) ne D then
	    "Solution does not have correct norm";
    end if;
    return quadratic_form, sol;
end function;

RunMeOnce := function(p, deg)
    B<i,j,k> := QuaternionAlgebra(p);
    O0 := QuaternionOrder( [ 1, i, 1/2*j + 1/2*i, 1/2 + 1/2*k ] );
    isoDegree := 2;
    quadratic_form, sol := Random_NormEqwSol(B, O0, deg);
    return quadratic_form, sol;
end function;

RunMeNTimes := function(p, deg, n)
    forms := [];
    sols := [];
    for i in [1..n] do
	form, sol := RunMeOnce(p, deg);
	forms[i] := form;
	sols[i] := sol;
    end for;
    return forms, sols;
end function;


