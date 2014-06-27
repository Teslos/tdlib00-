<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=pd format=gnuplot>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=Tm_s10></algorithm>
        <output name=x2> x(L,Se) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=s10_L_s32></algorithm>
        <output name=x2> x(L,Se) </output>
      </compute>
      <var name=x2 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=s10_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,Se) </input>
      <output> x(L,Se) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s10_L_s32></algorithm>
        <output name=x2> x(L,Se) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=s32_L_s11></algorithm>
        <output name=x2> x(L,Se) </output>
      </compute>
      <var name=x2 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=s32_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,Se) </input>
      <output> x(L,Se) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output ChangeXY>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s32_L_s11></algorithm>
        <output name=x2l> x(L,Se) </output>
        <output name=x2r> x(s11,Se) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=s11_L_s23></algorithm>
        <output name=T> T </output>
      </compute>
      <var name=T operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T + 5 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=s11_L_1></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2l> x(L,Se) </input>
      <input once name=x2r> x(s11,Se) </input>
      <output> T </output>
      <output noname> x(L,Se) </output>
      <output noname> x(s11,Se) </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s11_L_s23></algorithm>
        <output name=x2> x(L,Se) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=s23_L1_L2></algorithm>
        <output name=x2> x(L_1,2) </output>
      </compute>
      <var name=x2 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=s23_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,Se) </input>
      <output> x(L,Se) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output ChangeXY>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s23_L1_L2></algorithm>
        <output name=x2l> x(L_1,2) </output>
        <output name=x2r> x(L_2,2) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <var name=T value=1050 operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T + 5 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=L1_L2></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2l> x(L_1,2) </input>
      <input once name=x2r> x(L_2,2) </input>
      <output> T </output>
      <output noname> x(L_1,2) </output>
      <output noname> x(L_2,2) </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s23_L1_L2></algorithm>
        <output name=x2> x(L_2,2) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=s23_L_s01></algorithm>
        <output name=x2> x(L,Se) </output>
      </compute>
      <var name=x2 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=s23_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,Se) </input>
      <output> x(L,Se) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s23_L_s01></algorithm>
        <output name=x2> x(L,Se) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=s01_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,Se) </input>
      <output> x(L,Se) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s32_L_s11></algorithm>
        <output name=x2> x(s11,Se) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <var name=T value=400 operation=LE></var>
    </finish>
    <step>
      <var name=T><convert> T - 10 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=s11_s32></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2> x(s11,Se) </input>
      <output> x(s11,Se) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s11_L_s23></algorithm>
        <output name=x2> x(s11,Se) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <var name=T value=400 operation=LE></var>
    </finish>
    <step>
      <var name=T><convert> T - 15 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=s11_s23></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2> x(s11,Se) </input>
      <output> x(s11,Se) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s10_L_s32></algorithm>
        <output> T </output>
      </compute>
      <var name=x2 value=0></var>
    </start>
    <finish>
      <var name=x2 value=0.4 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 + 1 </convert></var>
    </step>
    <compute> 
      <PassThrough class=algorithm id=p5></PassThrough>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s32_L_s11></algorithm>
        <output name=x2> x(L,Se) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=s32_L_s11></algorithm>
        <output name=x2> x(s11,Se) </output>
      </compute>
    </finish>
    <step>
      <var name=x2><convert> x2 + 1 </convert></var>
    </step>
    <compute> 
      <PassThrough class=algorithm id=p6></PassThrough>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s11_L_s23></algorithm>
        <output name=x2> x(L,Se) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <var name=x2 value=0.6 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 + 1 </convert></var>
    </step>
    <compute> 
      <PassThrough class=algorithm id=p7></PassThrough>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s23_L1_L2></algorithm>
        <output> T </output>
      </compute>
      <var name=x2 value=0.6></var>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=s23_L1_L2></algorithm>
        <output name=x2> x(L_2,2) </output>
      </compute>
    </finish>
    <step>
      <var name=x2><convert> x2 + 1 </convert></var>
    </step>
    <compute> 
      <PassThrough class=algorithm id=p1></PassThrough>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=s23_L_s01></algorithm>
        <output> T </output>
      </compute>
      <var name=x2 value=0.6></var>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 + 1 </convert></var>
    </step>
    <compute> 
      <PassThrough class=algorithm id=p2></PassThrough>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <var name=T value=400></var>
      <var name=x2 value=0.4></var>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=s32_L_s11></algorithm>
        <output> T </output>
      </compute>
    </finish>
    <step>
      <var name=T><convert> T + 2000 </convert></var>
    </step>
    <compute> 
      <PassThrough class=algorithm id=p3></PassThrough>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <var name=T value=400></var>
      <var name=x2 value=0.6></var>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=Tm_s23></algorithm>
        <output> T </output>
      </compute>
    </finish>
    <step>
      <var name=T><convert> T + 2000 </convert></var>
    </step>
    <compute> 
      <PassThrough class=algorithm id=p4></PassThrough>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
</OutputFile>

