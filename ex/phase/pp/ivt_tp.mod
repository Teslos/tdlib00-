<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=test>
  <species>
    <IVT_Tp class=func_Tp>
      (
      <coef> 10 </coef> +
      <coef> 5 </coef> *T+
      <coef> -2 </coef> *T*log(T)+
      <coef> -1e-5 </coef> *T^2+
      <coef> -2e-8 </coef> *T^3+
      <coef> -3e-12 </coef> *T^4+
      <coef> 2000 </coef> /T)
    </IVT_Tp> 
  </species>
</PointPhase>

