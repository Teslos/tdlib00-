<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=t format=file>
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=1000 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T+300 </convert>
      </var>
    </step>
    <compute>
      <PhaseProperty class=algorithm id=test
        DependentMoleFraction=first>
        <phase class=phase IDREF=test></phase>
      </PhaseProperty> 
      <input value=1> p </input>
      <input name=T> T </input>
      <output> T </output>
      <output> G </output>
      <output> G_ref </output>
      <output> G_mix </output>
      <output> G_excess </output>
      <output> G_ideal </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=1000 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T+300 </convert>
      </var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=test></algorithm>
      <input value=1> p </input>
      <input name=T> T </input>
      <output> T </output>
      <output> H </output>
      <output> H_ref </output>
      <output> H_mix </output>
      <output> H_excess </output>
      <output> H_ideal </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=1000 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T+300 </convert>
      </var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=test></algorithm>
      <input value=1> p </input>
      <input name=T> T </input>
      <output> T </output>
      <output> S </output>
      <output> S_ref </output>
      <output> S_mix </output>
      <output> S_excess </output>
      <output> S_ideal </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

