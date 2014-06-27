<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<SimpleSolution class=phase id=test>
  <components> _A _B _C </components>
  <HochArpshofen class=FuncTpx formalism=NoChange>
    +x(_A)*x(_B)/v(_A)/v(_B)*(
      <func_x> +(v(_B)-v(_B)^2)*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> -9000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +(v(_B)-v(_B)^3)*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> -1000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </HochArpshofen> 
</SimpleSolution>
<SimpleSolution class=phase id=testM>
  <components> _A _B _C </components>
  <HochArpshofen class=FuncTpx formalism=Muggianu>
    +x(_A)*x(_B)/v(_A)/v(_B)*(
      <func_x> +(v(_B)-v(_B)^2)*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> -9000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +(v(_B)-v(_B)^3)*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> -1000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </HochArpshofen> 
</SimpleSolution>
<SimpleSolution class=phase id=testK>
  <components> _A _B _C </components>
  <HochArpshofen class=FuncTpx formalism=Kohler>
    +x(_A)*x(_B)/v(_A)/v(_B)*(
      <func_x> +(v(_B)-v(_B)^2)*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> -9000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +(v(_B)-v(_B)^3)*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> -1000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </HochArpshofen> 
</SimpleSolution>

