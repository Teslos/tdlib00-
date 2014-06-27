<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=frml format=file>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0></var>
      <var name=x3 value=0></var>
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
      <input name=x3> x(_C) </input>
      <output> x(_A) </output>
      <output> x(_B) </output>
      <output name=G_NoChange> G </output>
    </compute>
    <compute>
      <algorithm IDREF=testM></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(_B) </input>
      <input name=x3> x(_C) </input>
      <output name=G_Muggianu> G </output>
    </compute>
    <compute>
      <algorithm IDREF=testK></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(_B) </input>
      <input name=x3> x(_C) </input>
      <output name=G_Kohler> G </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0></var>
      <var name=x3 value=0.1></var>
    </start>
    <finish>
      <var name=x2 value=0.9 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.2 </convert>
      </var>
    </step>
    <compute>
      <algorithm IDREF=test></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(_B) </input>
      <input name=x3> x(_C) </input>
      <output> x(_A) </output>
      <output> x(_B) </output>
      <output name=G_NoChange> G </output>
    </compute>
    <compute>
      <algorithm IDREF=testM></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(_B) </input>
      <input name=x3> x(_C) </input>
      <output name=G_Muggianu> G </output>
    </compute>
    <compute>
      <algorithm IDREF=testK></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(_B) </input>
      <input name=x3> x(_C) </input>
      <output name=G_Kohler> G </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

