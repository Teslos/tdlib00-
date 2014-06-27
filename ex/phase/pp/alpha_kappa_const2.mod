<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=test>
  <species>
    <alpha_kappa_const2 class=func_Tp>
      (
      <coef> 100 </coef> /kappa*exp(
      <coef> 1e-4 </coef> *T)*(1.-exp(
      <coef> 1e-5 </coef> *(1-p/po)))
    </alpha_kappa_const2> 
  </species>
</PointPhase>

