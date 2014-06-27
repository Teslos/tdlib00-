<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PhaseProperty class=algorithm id=test
  DependentMoleFraction=first>
  <phase class=phase IDREF=test></phase>
</PhaseProperty> 

<OutputFile ext=x1 format=file>
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
      <algorithm IDREF=test></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(_B) </input>
      <output> x(_B) </output>
      <output> G </output>
      <output> H </output>
      <output> S </output>
      <output> Cp </output>
    </compute>
  </ComputeOutput> 
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
      <algorithm IDREF=test></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(_B) </input>
      <output> x(_B) </output>
      <output> G(all) </output>
      <output> H(all) </output>
    </compute>
  </ComputeOutput> 
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
      <algorithm IDREF=test></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(_B) </input>
      <output> x(_B) </output>
      <output> G_ref </output>
      <output> G_mix </output>
      <output> G_ideal </output>
      <output> G_excess </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

