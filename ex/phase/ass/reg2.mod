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
    <components> Hg Te Hg0.5Te0.5 Hg0.2Te0.8 </components>
    <Reference class=FuncTpx>
      <func_x> +x(Hg)* </func_x>
      <species> Hg
      </species>
      <func_x> +x(Te)* </func_x>
      <species> Te
      </species>
      <func_x> +x(Hg0.5Te0.5)* </func_x>
      <species> Hg0.5Te0.5
        <Cp_zero class=func_Tp>
          (
          <coef> -10000 </coef> +
          <coef> 0 </coef> *T)
        </Cp_zero> 
      </species>
      <func_x> +x(Hg0.2Te0.8)* </func_x>
      <species> Hg0.2Te0.8
        <Cp_zero class=func_Tp>
          (
          <coef> -30000 </coef> +
          <coef> 0 </coef> *T)
        </Cp_zero> 
      </species>
    </Reference> 
    <IdealMixing class=FuncTpx>
    </IdealMixing>
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Hg)*x(Hg0.5Te0.5)*(
        <func_x> + </func_x> 
        <Cp_zero class=func_Tp>
          (
          <coef> -5400 </coef> +
          <coef> 1.2 </coef> *T)
        </Cp_zero> 
      )
    </RedlichKister> 
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Hg)*x(Hg0.2Te0.8)*(
        <func_x> + </func_x> 
        <Cp_zero class=func_Tp>
          (
          <coef> -3300 </coef> +
          <coef> 2.2 </coef> *T)
        </Cp_zero> 
      )
    </RedlichKister> 
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Te)*x(Hg0.5Te0.5)*(
        <func_x> + </func_x> 
        <Cp_zero class=func_Tp>
          (
          <coef> +1100 </coef> +
          <coef> 0.5 </coef> *T)
        </Cp_zero> 
      )
    </RedlichKister> 
    <RedlichKister class=FuncTpx formalism=Muggianu>
      +x(Te)*x(Hg0.2Te0.8)*(
        <func_x> + </func_x> 
        <Cp_zero class=func_Tp>
          (
          <coef> +3300 </coef> +
          <coef> 0.8 </coef> *T)
        </Cp_zero> 
      )
    </RedlichKister> 
  </SimpleSolution> 
</AssociatedSolution> 


