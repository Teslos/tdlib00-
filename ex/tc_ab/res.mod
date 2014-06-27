<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<residual id=BCC_L_A2B_T
  yname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BCC_L_A2B></algorithm>
    <output> T </output>
  </compute>
</residual>
<residual id=BCC_L_A2B_L
  yname=x2l
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BCC_L_A2B></algorithm>
    <output> x(L,B) </output>
  </compute>
</residual>
<residual id=BCC_L_A2B_B
  yname=x2bcc
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BCC_L_A2B></algorithm>
    <output> x(BCC,B) </output>
  </compute>
</residual>

<residual id=A2B_L
  yname=T
  xname=x2
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=A2B_L></algorithm>
    <input name=x2> x(L,B) </input>
    <input name=T> T </input>
    <output> T </output>
  </compute>
</residual>

<residual id=A2B_L_BCC_T
  yname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=A2B_L_BCC></algorithm>
    <output> T </output>
  </compute>
</residual>
<residual id=A2B_L_BCC_L
  yname=x2l
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=A2B_L_BCC></algorithm>
    <output> x(L,B) </output>
  </compute>
</residual>
<residual id=A2B_L_BCC_B
  yname=x2bcc
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=A2B_L_BCC></algorithm>
    <output> x(BCC,B) </output>
  </compute>
</residual>

<residual id=L_BCC_FCC_T
  yname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
    <output> T </output>
  </compute>
</residual>
<residual id=L_BCC_FCC_L
  yname=x2l
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
    <output> x(L,B) </output>
  </compute>
</residual>
<residual id=L_BCC_FCC_B
  yname=x2bcc
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
    <output> x(BCC,B) </output>
  </compute>
</residual>
<residual id=L_BCC_FCC_F
  yname=x2fcc
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
    <output> x(FCC,B) </output>
  </compute>
</residual>

<residual id=BCC_A2B_BCC_T
  yname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
    <output> T </output>
  </compute>
</residual>
<residual id=BCC_A2B_BCC_B1
  yname=x2_1
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
    <output> x(BCC_1,B) </output>
  </compute>
</residual>
<residual id=BCC_A2B_BCC_B2
  yname=x2_2
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
    <output> x(BCC_2,B) </output>
  </compute>
</residual>

<residual id=FCC_L_L1
  yname=x2l
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=FCC_L></algorithm>
    <input name=x2l> x(L,B) </input>
    <input value=0.9> x(FCC,B) </input>
    <input name=T> T </input>
    <output> x(L,B) </output>
  </compute>
</residual>
<residual id=FCC_L_L2
  yname=x2l
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=FCC_L></algorithm>
    <input name=x2l> x(L,B) </input>
    <input name=x2fcc> x(FCC,B) </input>
    <input name=T> T </input>
    <output> x(L,B) </output>
  </compute>
</residual>
<residual id=FCC_L_F
  yname=x2fcc
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=FCC_L></algorithm>
    <input name=x2l> x(L,B) </input>
    <input name=x2fcc> x(FCC,B) </input>
    <input name=T> T </input>
    <output> x(FCC,B) </output>
  </compute>
</residual>

<residual id=BCC_L_L
  yname=x2l
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BCC_L_r></algorithm>
    <input name=x2l> x(L,B) </input>
    <input name=x2bcc> x(BCC,B) </input>
    <input name=T> T </input>
    <output> x(L,B) </output>
  </compute>
</residual>
<residual id=BCC_L_B
  yname=x2bcc
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=BCC_L_r></algorithm>
    <input name=x2l> x(L,B) </input>
    <input name=x2bcc> x(BCC,B) </input>
    <input name=T> T </input>
    <output> x(BCC,B) </output>
  </compute>
</residual>

<residual id=act
  yname=aB
  xname=x2
  ScaleOfX=1>
  <compute>
    <reaction class=algorithm id=r1>
      <compute>
        <algorithm class=algorithm IDREF=L></algorithm>
        <input once name=T> T </input>
        <input name=x2> x(B) </input>
        <output name=mu> G(B) </output>
        <output> T </output>
      </compute>
      <compute>
        <algorithm class=algorithm IDREF=Bfcc></algorithm>
        <input once name=T> T </input>
        <output name=muo> G </output>
      </compute>
    </reaction>
    <input name=x2> x2 </input>
    <input once name=T> T </input>
    <output> all </output>
    <convert> exp((mu-muo)/R/T) </convert>
  </compute>
</residual>

<residual id=Hmix
  yname=Hmix
  xname=x2
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=L></algorithm>
    <input name=x2> x(B) </input>
    <input once name=T> T </input>
    <output> H_mix </output>
  </compute>
</residual>

<residual id=Hmelt
  yname=H
  xname=T
  ScaleOfX=1>
  <compute>
    <algorithm class=algorithm IDREF=A2B_L></algorithm>
    <input name=x2> x(L,B) </input>
    <input name=T> T </input>
    <output name=A2B> H(A2B) </output>
    <output name=L> H(L) </output>
    <convert> L-A2B/3 </convert>
  </compute>
</residual>

