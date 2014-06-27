<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=test>
  <species>
    <alpha_const class=func_Tp>
      (
      <coef> 100 </coef> *exp(
      <coef> 1e-4 </coef> *T)*(p/po-1)
    </alpha_const> 
  </species>
</PointPhase>

