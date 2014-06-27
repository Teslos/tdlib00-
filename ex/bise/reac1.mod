<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=reac1 format=file>
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
      <reaction class=algorithm id=reac1>
        <compute>
          <algorithm class=algorithm IDREF=s23></algorithm>
          <input name=T> T </input>
          <output> T </output>
          <output name=s23> G </output>
        </compute>
        <compute>
          <algorithm class=algorithm IDREF=s10></algorithm>
          <input name=T> T </input>
          <output name=s10> G </output>
        </compute>
        <compute>
          <algorithm class=algorithm IDREF=s01></algorithm>
          <input name=T> T </input>
          <output name=s01> G </output>
        </compute>
      </reaction> 
      <input name=T> T </input>
      <output> all </output>
      <convert name=T> T </convert>
      <convert name=EMF> -(5*s23 - 2*s10 - 3*s01)/(6*96.520) </convert>
    </compute>
  </ComputeOutput>
</OutputFile>
