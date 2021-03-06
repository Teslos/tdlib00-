<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=test>
  <species>
    <Cp_BB4 class=func_Tp>
      (
      <coef> 10 </coef> +
      <coef> 5 </coef> *T+
      <coef> -2 </coef> *T*log(T)+
      <coef> -0.5 </coef> *sqrt(T)+
      <coef> 2500 </coef> /T+
      <coef> 500000 </coef> /T^2)
    </Cp_BB4> 
  </species>
</PointPhase>

