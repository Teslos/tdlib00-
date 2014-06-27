<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<AssociatedSolution class=phase id=L
  equilibrated=0 debug=0>
  <SimpleSolution class=phase>
    <components> Hg Te </components>
    <Reference class=FuncTpx>
      <func_x> +x(Hg)* </func_x>
      <species> Hg
      </species>
      <func_x> +x(Te)* </func_x>
      <species> Te
      </species>
    </Reference> 
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
          <coef> -42000 </coef> +
          <coef> 2.5 </coef> *T)
        </Cp_zero> 
      </species>
    </Reference> 
    <IdealMixing class=FuncTpx>
      +R*T*x(Hg)*log(x(Hg))
      +R*T*x(Te)*log(x(Te))
      +R*T*x(HgTe)*log(x(HgTe))
    </IdealMixing> 
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Hg)*x(HgTe)*(
        <func_x> + </func_x> 
        <Cp_zero class=func_Tp>
          (
          <coef> -5400 </coef> +
          <coef> 1.2 </coef> *T)
        </Cp_zero> 
      )
    </RedlichKister> 
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Te)*x(HgTe)*(
        <func_x> + </func_x> 
        <Cp_zero class=func_Tp>
          (
          <coef> +1100 </coef> +
          <coef> 0.5 </coef> *T)
        </Cp_zero> 
      )
    </RedlichKister> 
  </SimpleSolution> 
</AssociatedSolution> 


