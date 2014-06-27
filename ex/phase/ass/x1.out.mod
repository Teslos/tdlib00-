<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=x1>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0.001></var>
    </start>
    <finish>
      <var name=x2 value=0.999 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.2 </convert>
      </var>
    </step>
    <compute>
      <algorithm IDREF=twosol></algorithm>
      <input once value=1000> T </input>
      <input name=x2> x(2) </input>
      <output> x(2) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

