<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<SimpleSolution class=phase id=test>
  <components> _A _B _C </components>
  <Polynomial class=FuncTpx formalism=NoChange>
    +x(_A)*x(_B)*(
      <func_x> +v(1)^2*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 10000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +v(1)*v(2)*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 19000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +v(2)^2*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 9100 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </Polynomial> 
</SimpleSolution>
<SimpleSolution class=phase id=testM>
  <components> _A _B _C </components>
  <Polynomial class=FuncTpx formalism=Muggianu>
    +x(_A)*x(_B)*(
      <func_x> +v(1)^2*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 10000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +v(1)*v(2)*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 19000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +v(2)^2*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 9100 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </Polynomial> 
</SimpleSolution>
<SimpleSolution class=phase id=testK>
  <components> _A _B _C </components>
  <Polynomial class=FuncTpx formalism=Kohler>
    +x(_A)*x(_B)*(
      <func_x> +v(1)^2*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 10000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +v(1)*v(2)*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 19000 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
      <func_x> +v(2)^2*</func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef> 9100 </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </Polynomial> 
</SimpleSolution>

