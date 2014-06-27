<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=small format=file>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=-0.004></var>
    </start>
    <finish>
      <var name=x2 value=0 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.001 </convert>
      </var>
    </step>
    <compute>
      <algorithm IDREF=test></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(_B) </input>
      <output> x(_B) </output>
      <output> G </output>
      <output> G(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <var name=lgx value=-15></var>
    </start>
    <finish>
      <var name=lgx value=-5 operation=GE></var>
    </finish>
    <step>
      <var name=lgx><convert> lgx+1 </convert>
      </var>
    </step>
    <compute>
      <algorithm IDREF=test></algorithm>
      <input value=1000> T </input>
      <input><convert>10^lgx</convert> x(_B) </input>
      <output> x(_B) </output>
      <output> G </output>
      <output> G(all) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

