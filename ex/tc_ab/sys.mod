<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<elements>
  A    20
  B    50
</elements>
  
<PointPhase class=phase id=A2B>
  <species> A2B
    <Cp_const class=func_Tp>
      (
      <coef id=V1 unknown=1 scale=50000> 20276.46008 </coef> +
      <coef id=V2 unknown=1 scale=50> -29.54712198 </coef> *T+
      <coef id=V3 unknown=0 scale=50> 0 </coef> *T*log(T))
    </Cp_const> 
  </species>
</PointPhase> 
<SimpleSolution class=phase id=BCC>
  <components> A B </components>
  <Reference class=FuncTpx>
    <func_x> +x(A)* </func_x>
    <species id=A_bcc> A
    </species>
    <func_x> +x(B)* </func_x>
    <species id=B_bcc> B
    </species>
  </Reference> 
  <IdealMixing class=FuncTpx>
    +R*T*x(A)*log(x(A))
    +R*T*x(B)*log(x(B))
  </IdealMixing> 
  <RedlichKister class=FuncTpx formalism=NoChange>
    +x(A)*x(B)*(
      <func_x> + </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef id=V15 unknown=1 scale=50000> 25199.38702 </coef> +
        <coef id=V16 unknown=1 scale=50> -9.697850716 </coef> *T)
      </Cp_zero> 
      <func_x> +(x(A)-x(B))* </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef id=V17 unknown=1 scale=50000> 3906.021158 </coef> +
        <coef id=V18 unknown=0 scale=50> 0 </coef> *T)
      </Cp_zero> 
    )
  </RedlichKister> 
</SimpleSolution> 
<SimpleSolution class=phase id=FCC>
  <components> A B </components>
  <Reference class=FuncTpx>
    <func_x> +x(A)* </func_x>
    <species> A
      <Cp_zero class=func_Tp>
        (
        <coef> 408 </coef> +
        <coef> 0 </coef> *T)
      </Cp_zero> 
    </species>
    <func_x> +x(B)* </func_x>
    <species id=B_fcc> B
      <Cp_zero class=func_Tp>
        (
        <coef> 3300 </coef> +
        <coef> -3 </coef> *T)
      </Cp_zero> 
    </species>
  </Reference> 
  <IdealMixing class=FuncTpx>
    +R*T*x(A)*log(x(A))
    +R*T*x(B)*log(x(B))
  </IdealMixing> 
  <RedlichKister class=FuncTpx formalism=NoChange>
    +x(A)*x(B)*(
      <func_x> + </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef id=V19 unknown=1 scale=50000> 25000 </coef> +
        <coef id=V20 unknown=1 scale=50> -9.7 </coef> *T)
      </Cp_zero> 
      <func_x> +(x(A)-x(B))* </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef id=V21 unknown=0 scale=50000> 0 </coef> +
        <coef id=V22 unknown=0 scale=50> 0 </coef> *T)
      </Cp_zero> 
    )
  </RedlichKister> 
</SimpleSolution> 
<SimpleSolution class=phase id=L>
  <components> A B </components>
  <Reference class=FuncTpx>
    <func_x> +x(A)* </func_x>
    <species id=A_l> A
      <Cp_zero class=func_Tp>
        (
        <coef> 14000 </coef> +
        <coef> -10 </coef> *T)
      </Cp_zero> 
    </species>
    <func_x> +x(B)* </func_x>
    <species id=B_l> B
      <Cp_zero class=func_Tp>
        (
        <coef> 18000 </coef> +
        <coef> -12 </coef> *T)
      </Cp_zero> 
    </species>
  </Reference> 
  <IdealMixing class=FuncTpx>
    +R*T*x(A)*log(x(A))
    +R*T*x(B)*log(x(B))
  </IdealMixing> 
  <RedlichKister class=FuncTpx formalism=NoChange>
    +x(A)*x(B)*(
      <func_x> + </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef id=V11 unknown=1 scale=50000> -21813.80981 </coef> +
        <coef id=V12 unknown=1 scale=50> 14.97997819 </coef> *T)
      </Cp_zero> 
      <func_x> +(x(A)-x(B))* </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef id=V13 unknown=0 scale=50000> 0 </coef> +
        <coef id=V14 unknown=0 scale=50> 0 </coef> *T)
      </Cp_zero> 
    )
  </RedlichKister> 
</SimpleSolution> 
