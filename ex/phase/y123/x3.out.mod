<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=x3>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0></var>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.2 </convert>
      </var>
    </step>
    <compute>
      <algorithm IDREF=twosol></algorithm>
      <input once value=1000> T </input>
      <input name=x2> x(YBa2Cu3O7) </input>
      <output> x(YBa2Cu3O7) </output>
      <output> S(all) </output>
      <output> V(all) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

