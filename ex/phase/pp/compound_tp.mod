<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=test>
  <species> 
    <compound_Tp class=func_Tp>
      <Cp_BB2 class=func_Tp>
        (
        <coef> 10 </coef> +
        <coef> 5 </coef> *T+
        <coef> -2 </coef> *T*log(T)+
        <coef> -0.4 </coef> *sqrt(T))
      </Cp_BB2> 
      <limit_Tp> T = {200, 600} p = {0, inf} </limit_Tp>
      <ideal_gas class=func_Tp>
        (R*T*log(p))
      </ideal_gas> 
      <limit_Tp> T = {600, 1500} p = {0, inf} </limit_Tp>
    </compound_Tp>   </species>
</PointPhase>

