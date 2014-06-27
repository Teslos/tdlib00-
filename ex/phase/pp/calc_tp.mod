<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PointPhase class=phase id=test>
  <species>
    <calc_Tp class=func_Tp>
      (
       10 + 5*T - 2*T*log(T) - 0.4*sqrt(T)) + R*T*log(p)
      )
    </calc_Tp> 
  </species>
</PointPhase>

