<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<SimpleSolution class=phase id=test>
  <components> _A _B </components>
  <IdealMixing>
    +R*T*<coef>2</coef>*x(_A)*log(x(_A))
    +R*T*<coef>3</coef>*x(_B)*log(x(_B))
  </IdealMixing>
</SimpleSolution>

