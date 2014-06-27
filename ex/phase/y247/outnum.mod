<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=num format=file>
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
      <algorithm IDREF=Y247></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(Y2Ba4Cu7O15) </input>
      <output> x(Y2Ba4Cu7O15) </output>
      <output> H </output>
    </compute>
    <compute>
      <algorithm IDREF=Y247_num></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(Y2Ba4Cu7O15) </input>
      <output name=H_num> H </output>
    </compute>
    <compute>
      <algorithm IDREF=Y247></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(Y2Ba4Cu7O15) </input>
      <output> Cp </output>
    </compute>
    <compute>
      <algorithm IDREF=Y247_num></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(Y2Ba4Cu7O15) </input>
      <output name=Cp_num> Cp </output>
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
      <algorithm IDREF=Y247></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(Y2Ba4Cu7O15) </input>
      <output> x(Y2Ba4Cu7O15) </output>
      <output name=G1> G(Y2Ba4Cu7O14) </output>
    </compute>
    <compute>
      <algorithm IDREF=Y247_num></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(Y2Ba4Cu7O15) </input>
      <output name=G1_num> G(Y2Ba4Cu7O14) </output>
    </compute>
    <compute>
      <algorithm IDREF=Y247></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(Y2Ba4Cu7O15) </input>
      <output name=G2> G(Y2Ba4Cu7O15) </output>
    </compute>
    <compute>
      <algorithm IDREF=Y247_num></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(Y2Ba4Cu7O15) </input>
      <output name=G2_num> G(Y2Ba4Cu7O15) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

