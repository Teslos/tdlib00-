<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<SimpleSolution class=phase id=test>
  <components> _A _B _C </components>
  <RedlichKister class=FuncTpx formalism=NoChange>
    +x(_A)*x(_B)*(
      <func_x> + </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> -10500 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +(x(_A)-x(_B))* </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 500 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </RedlichKister> 
</SimpleSolution>
<SimpleSolution class=phase id=testM>
  <components> _A _B _C </components>
  <RedlichKister class=FuncTpx formalism=Muggianu>
    +x(_A)*x(_B)*(
      <func_x> + </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> -10500 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +(x(_A)-x(_B))* </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 500 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </RedlichKister> 
</SimpleSolution>
<SimpleSolution class=phase id=testK>
  <components> _A _B _C </components>
  <RedlichKister class=FuncTpx formalism=Kohler>
    +x(_A)*x(_B)*(
      <func_x> + </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> -10500 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +(x(_A)-x(_B))* </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 500 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </RedlichKister> 
</SimpleSolution>

