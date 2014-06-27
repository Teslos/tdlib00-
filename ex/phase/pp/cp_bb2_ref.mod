<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=test>
  <species>
    <Cp_BB2_Tref class=func_Tp
      To=300>
      (
      <coef id=s11_lHo unknown=1> 27513.8 </coef> -
      <coef id=s11_lSb unknown=1> 14.8122 </coef> *T+
      <coef id=s11_Cpa unknown=1> 4.01738 </coef> *(T-To-T*log(T/To))+4*
      <coef id=s11_Cpb unknown=1> -34.1601 </coef> *(sqrt(T)-sqrt(To)))
    </Cp_BB2_Tref> 
  </species>
</PointPhase>

