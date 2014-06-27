<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=reac2 format=file>
  <ComputeOutput>
    <start>
      <var name=T value=300></var>
    </start>
    <finish>
      <var name=T value=600 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T+50 </convert></var>
    </step>
    <compute>
      <reaction class=algorithm id=reac2>
        <compute>
          <algorithm class=algorithm IDREF=s23></algorithm>
          <input name=T> T </input>
          <output> T </output>
          <output name=HT> H </output>
        </compute>
        <compute>
          <algorithm class=algorithm IDREF=s23></algorithm>
          <input value=298.15> T </input>
          <output name=H298> H </output>
        </compute>
      </reaction> 
      <input name=T> T </input>
      <output> all </output>
      <convert name=T> T </convert>
      <convert name=HT_H298> (HT-H298)/1000 </convert>
    </compute>
  </ComputeOutput>
</OutputFile>
