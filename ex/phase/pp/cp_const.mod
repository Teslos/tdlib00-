<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=test>
  <species>
    <Cp_const class=func_Tp>
      (
      <coef> 10 </coef> +
      <coef> 5 </coef> *T+
      <coef> -1.5 </coef> *T*log(T))
    </Cp_const> 
  </species>
</PointPhase>

