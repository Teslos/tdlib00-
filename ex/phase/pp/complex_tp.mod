<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=test>
  <species>
    <complex_Tp class=func_Tp>
      (
      <Cp_BB2 class=func_Tp>
        (
        <coef> 10 </coef> +
        <coef> 5 </coef> *T+
        <coef> -2 </coef> *T*log(T)+
        <coef> -0.4 </coef> *sqrt(T))
      </Cp_BB2> 
      +
      <ideal_gas class=func_Tp>
        (R*T*log(p))
      </ideal_gas> 
      )
    </complex_Tp> 
  </species>
</PointPhase>

