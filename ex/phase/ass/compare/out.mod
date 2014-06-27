<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=int format=gnuplot>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0.001></var>
    </start>
    <finish>
      <var name=x2 value=0.999 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.0025 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=Leq></algorithm>
      <input once value=1000> T </input>
      <input name=x2> x(2) </input>
      <output> x(2) </output>
      <output name=Hg_eq> internal(x(Hg)) </output>
      <output name=Te_eq> internal(x(Te)) </output>
      <output name=HgTe_eq> internal(x(HgTe)) </output>
    </compute>
    <compute>
      <algorithm IDREF=Lnew></algorithm>
      <input once value=1000> T </input>
      <input name=x2> x(2) </input>
      <output name=Hg_new> internal(x(Hg)) </output>
      <output name=Te_new> internal(x(Te)) </output>
      <output name=HgTe_new> internal(x(HgTe)) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>
<OutputFile ext=g format=gnuplot>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0.01></var>
    </start>
    <finish>
      <var name=x2 value=0.99 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.0025 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=Leq></algorithm>
      <input once value=1000> T </input>
      <input name=x2> x(2) </input>
      <output> x(2) </output>
      <output name=G_eq> G </output>
      <output name=GHg_eq> G(Hg) </output>
      <output name=GTe_eq> G(Te) </output>
    </compute>
    <compute>
      <algorithm IDREF=Lnew></algorithm>
      <input once value=1000> T </input>
      <input name=x2> x(2) </input>
      <output name=G_new> G </output>
      <output name=GHg_new> G(Hg) </output>
      <output name=GTe_new> G(Te) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

