<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<SimpleSolution class=phase id=test>
  <components> _A _B </components>
  <Reference class=FuncTpx>
    <func_x> +x(_A)* </func_x>
    <species>
      <Cp_zero class=func_Tp>
        (
        <coef> 1000 </coef> +
        <coef> -0.5 </coef> *T)
      </Cp_zero> 
    </species>
    <func_x> +x(_B)* </func_x>
    <species>
      <Cp_zero class=func_Tp>
        (
        <coef> 10000 </coef> +
        <coef> -5 </coef> *T)
      </Cp_zero> 
    </species>
  </Reference>
</SimpleSolution>

