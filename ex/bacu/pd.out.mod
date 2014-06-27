<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=pd format=gnuplot>
  <SeriesOutput Residuals=0> L1 </SeriesOutput>
  <SeriesOutput Residuals=0> L2 </SeriesOutput>
  <SeriesOutput Residuals=0> L3 </SeriesOutput>
  <SeriesOutput Residuals=0> L4 </SeriesOutput>
  <SeriesOutput Residuals=0> L5 </SeriesOutput>
  <SeriesOutput Residuals=0> L6 </SeriesOutput>
  <SeriesOutput Residuals=0> L7 </SeriesOutput>
  <SeriesOutput Residuals=0> L8 </SeriesOutput>
  <SeriesOutput Residuals=0> N1 </SeriesOutput>
  <SeriesOutput Residuals=0> N2 </SeriesOutput>
  <SeriesOutput Residuals=0> N3 </SeriesOutput>
  <SeriesOutput Residuals=0> N4 </SeriesOutput>
  <SeriesOutput Residuals=0> N5 </SeriesOutput>
  <SeriesOutput Residuals=0> N6 </SeriesOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=TmBa></algorithm>
        <output name=x2> x(L,Cu) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=Ba_L_BaCu></algorithm>
        <output name=x2> x(L,Cu) </output>
      </compute>
      <var name=x2 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=Ba_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,Cu) </input>
      <output> x(L,Cu) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=Ba_L_BaCu></algorithm>
        <output name=x2> x(L,Cu) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=BaCu_L_BaCu13></algorithm>
        <output name=x2> x(L,Cu) </output>
      </compute>
      <var name=x2 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute> <algorithm class=algorithm IDREF=BaCu_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,Cu) </input>
      <output> x(L,Cu) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BaCu_L_BaCu13></algorithm>
        <output name=x2> x(L,Cu) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=BaCu13_L_Cu></algorithm>
        <output name=x2> x(L,Cu) </output>
      </compute>
      <var name=x2 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute> <algorithm class=algorithm IDREF=BaCu13_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,Cu) </input>
      <output> x(L,Cu) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BaCu13_L_Cu></algorithm>
        <output name=x2> x(L,Cu) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute> 
      <algorithm class=algorithm IDREF=Cu_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,Cu) </input>
      <output> x(L,Cu) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=Ba_L_BaCu></algorithm>
        <output> T </output>
      </compute>
      <var name=x2 value=0></var>
    </start>
    <finish>
      <var name=x2 value=0.5 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 + 1 </convert></var>
    </step>
    <compute>     
      <algorithm class=algorithm IDREF=pass></algorithm>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BaCu_L_BaCu13></algorithm>
        <output> T </output>
      </compute>
      <var name=x2 value=0.5></var>
    </start>
    <finish>
      <var name=x2 value=0.92857 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 + 1 </convert></var>
    </step>
    <compute> 
      <algorithm class=algorithm IDREF=pass></algorithm>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BaCu13_L_Cu></algorithm>
        <output> T </output>
        <output name=x2> x(L, 2) </output>
      </compute>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 + 1 </convert></var>
    </step>
    <compute> <algorithm class=algorithm IDREF=pass></algorithm>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
      <var name=x2 value=0.5></var>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=BaCu_L></algorithm>
        <input value=0.5> x(L,2) </input>
        <input value=-1> T </input>
        <output> T </output>
      </compute>
    </finish>
    <step>
      <var name=T><convert> T + 2000 </convert></var>
    </step>
    <compute> 
      <algorithm class=algorithm IDREF=pass></algorithm>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
      <var name=x2 value=0.92857></var>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=BaCu13_L_Cu></algorithm>
        <output> T </output>
      </compute>
    </finish>
    <step>
      <var name=T><convert> T + 2000 </convert></var>
    </step>
    <compute> 
      <algorithm class=algorithm IDREF=pass></algorithm>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
</OutputFile>

