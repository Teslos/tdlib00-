<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<species id=Y> Y
  <Cp_zero class=func_Tp>
    (
    <coef> 2 </coef> +
    <coef> -3 </coef> *T)
  </Cp_zero>
</species>
<species id=Cu> Cu
  <Cp_zero class=func_Tp>
    (
    <coef> 4 </coef> +
    <coef> -6 </coef> *T)
  </Cp_zero>
</species>

<PointPhase class=phase id=test>
  <species> CuY
    <Cp_zero class=func_Tp>
      (
      <coef> 10 </coef> +
      <coef> -5 </coef> *T)
    </Cp_zero>
    <ref_plane>
      1
      <species IDREF=Cu></species>
      1
      <species IDREF=Y></species>
    </ref_plane>
  </species>
</PointPhase>

