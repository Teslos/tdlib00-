<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<AssociatedSolution class=phase id=Leq
  equilibrated=1 debug=0>
  <SimpleSolution class=phase>
    <components> Hg Te </components>
  </SimpleSolution> 
  <SimpleSolution class=phase>
    <components> Hg Te HgTe </components>
    <Reference class=FuncTpx>
      <func_x> +x(Hg)* </func_x>
      <species> Hg
      </species>
      <func_x> +x(Te)* </func_x>
      <species> Te
      </species>
      <func_x> +x(HgTe)* </func_x>
      <species> HgTe
        <Cp_zero class=func_Tp>
          (
          <coef id=ass> -50000 </coef> +
          <coef> 0 </coef> *T)
        </Cp_zero> 
      </species>
    </Reference> 
    <IdealMixing class=FuncTpx>
      +R*T*x(Hg)*log(x(Hg))
      +R*T*x(Te)*log(x(Te))
      +R*T*x(HgTe)*log(x(HgTe))
    </IdealMixing> 
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Te)*x(HgTe)*(
        <func_x> + </func_x> 
        <Cp_zero class=func_Tp>
          (
          <coef id=mis> 20000 </coef> +
          <coef> 0 </coef> *T)
        </Cp_zero> 
      )
    </RedlichKister> 
  </SimpleSolution> 
</AssociatedSolution> 

<AssociatedSolution class=phase id=Lnew
  equilibrated=0 debug=0>
  <SimpleSolution class=phase>
    <components> Hg Te </components>
  </SimpleSolution> 
  <SimpleSolution class=phase>
    <components> Hg Te HgTe </components>
    <Reference class=FuncTpx>
      <func_x> +x(Hg)* </func_x>
      <species> Hg
      </species>
      <func_x> +x(Te)* </func_x>
      <species> Te
      </species>
      <func_x> +x(HgTe)* </func_x>
      <species> HgTe
        <Cp_zero class=func_Tp>
          (
          <coef IDREF=ass></coef> +
          <coef> 0 </coef> *T)
        </Cp_zero> 
      </species>
    </Reference> 
    <IdealMixing class=FuncTpx>
      +R*T*x(Hg)*log(x(Hg))
      +R*T*x(Te)*log(x(Te))
      +R*T*x(HgTe)*log(x(HgTe))
    </IdealMixing> 
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Te)*x(HgTe)*(
        <func_x> + </func_x> 
        <Cp_zero class=func_Tp>
          (
          <coef IDREF=mis></coef> +
          <coef> 0 </coef> *T)
        </Cp_zero> 
      )
    </RedlichKister> 
  </SimpleSolution> 
</AssociatedSolution> 

