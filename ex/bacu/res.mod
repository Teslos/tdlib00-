<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<residual id=HBaCu
  yname=HBaCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BaCu></algorithm>
    <input value=298> T </input>
    <output> H_mix </output>
  </compute>
</residual>
<residual id=HBa
  yname=HBa
  xname=xCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=L></algorithm>
    <input name=xCu> x(Cu) </input>
    <input name=T> T </input>
    <output> H_mix(Ba) </output>
  </compute>
</residual>
<residual id=HCu
  yname=HCu
  xname=xCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=L></algorithm>
    <input name=xCu> x(Cu) </input>
    <input name=T> T </input>
    <output> H_mix(Cu) </output>
  </compute>
</residual>

<residual id=Ba_L
  yname=T
  xname=xCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=Ba_L></algorithm>
    <input name=T> T </input>
    <input name=xCu> x(L,Cu) </input>
    <output> T </output>
  </compute>
</residual>
<residual id=BaCu_L
  yname=T
  xname=xCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BaCu_L></algorithm>
    <input name=T> T </input>
    <input name=xCu> x(L,Cu) </input>
    <output> T </output>
  </compute>
</residual>
<residual id=BaCu13_L
  yname=T
  xname=xCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BaCu13_L></algorithm>
    <input name=T> T </input>
    <input name=xCu> x(L,Cu) </input>
    <output> T </output>
  </compute>
</residual>
<residual id=Cu_L
  yname=T
  xname=xCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=Cu_L></algorithm>
    <input name=T> T </input>
    <input name=xCu> x(L,Cu) </input>
    <output> T </output>
  </compute>
</residual>

<residual id=Ba_L_BaCu
  yname=T
  xname=xCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=Ba_L_BaCu></algorithm>
    <output> T </output>
  </compute>
</residual>
<residual id=BaCu_L_BaCu13
  yname=T
  xname=xCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BaCu_L_BaCu13></algorithm>
    <output> T </output>
  </compute>
</residual>
<residual id=BaCu13_L_Cu
  yname=T
  xname=xCu
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BaCu13_L_Cu></algorithm>
    <output> T </output>
  </compute>
</residual>

