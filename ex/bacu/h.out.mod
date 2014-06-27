<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=h format=gnuplot>
  <SeriesOutput Residuals=0> H1 </SeriesOutput>
  <SeriesOutput Residuals=0> H2 </SeriesOutput>
  <SeriesOutput Residuals=0> H3 </SeriesOutput>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0.0></var>
    </start>
    <finish>
      <var name=x2 value=1.0 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 + 0.01 </convert>
      </var>
    </step>
  <compute>
    <algorithm class=algorithm IDREF=L></algorithm>
    <input name=x2> x(Cu) </input>
    <input Value=1> p </input>
    <input Value=1400> T </input>
    <output> x(Cu) </output>
    <output> H_mix(Ba) </output>
    <output> H_mix(Cu) </output>
  </compute>
  </ComputeOutput>
</OutputFile>


