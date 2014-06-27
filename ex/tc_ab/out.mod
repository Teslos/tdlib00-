<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=non format=file>
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=TmA></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=TtrsB></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=TmB></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=BCC_L_A2B></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=A2B_L_BCC></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <compute>
      <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
      <output> state(all) </output>
      <output> fmin </output>
      <output> fmin(all) </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

<OutputFile ext=mis format=file>
  <SpinodalOutput phase=L direction=up></SpinodalOutput>
  <SpinodalOutput phase=L direction=down></SpinodalOutput>
  <SpinodalOutput phase=BCC direction=up></SpinodalOutput>
  <SpinodalOutput phase=BCC direction=down></SpinodalOutput>
  <SpinodalOutput phase=FCC direction=up></SpinodalOutput>
  <SpinodalOutput phase=FCC direction=down></SpinodalOutput>
</OutputFile>

<OutputFile ext=pd format=gnuplot>
  <SeriesOutput Residuals=0 ChangeXY> E1L </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> E1B </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> E2L </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> E2B </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> E3L </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> E3B </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> E3F </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> E4B1 </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> E4B2 </SeriesOutput>
  <SeriesOutput Residuals=0> L1 </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> L2 </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> L3L </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> L3F </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> L4L </SeriesOutput>
  <SeriesOutput Residuals=0 ChangeXY> L4B </SeriesOutput>
  <ComputeOutput class=output ChangeXY>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BCC_L_A2B></algorithm>
        <output name=x2l> x(L,B) </output>
        <output name=x2bcc> x(BCC,B) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=TmA></algorithm>
        <output name=T> T </output>
      </compute>
      <var name=T operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T + 5 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=BCC_L_l></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2l> x(L,B) </input>
      <input once name=x2bcc> x(BCC,B) </input>
      <output> T </output>
      <output noname> x(L,B) </output>
      <output noname> x(BCC,B) </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BCC_L_A2B></algorithm>
        <output name=x2> x(L,B) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=A2B_L_BCC></algorithm>
        <output name=x2> x(L,B) </output>
      </compute>
      <var name=x2 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=A2B_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input once name=T> T </input>
      <input name=x2> x(L,B) </input>
      <output> x(L,B) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output ChangeXY>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=A2B_L_BCC></algorithm>
        <output name=x2l> x(L,B) </output>
        <output name=x2bcc> x(BCC,B) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
        <output name=T> T </output>
      </compute>
      <var name=T operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T + 5 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=BCC_L_r></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2l> x(L,B) </input>
      <input once name=x2bcc> x(BCC,B) </input>
      <output> T </output>
      <output noname> x(L,B) </output>
      <output noname> x(BCC,B) </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output ChangeXY>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
        <output name=x2l> x(L,B) </output>
        <output name=x2fcc> x(FCC,B) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=TmB></algorithm>
        <output name=T> T </output>
      </compute>
      <var name=T operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T + 5 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=FCC_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2l> x(L,B) </input>
      <input once name=x2fcc> x(FCC,B) </input>
      <output> T </output>
      <output noname> x(L,B) </output>
      <output noname> x(FCC,B) </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output ChangeXY>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
        <output name=x2fcc> x(FCC,B) </output>
        <output name=x2bcc> x(BCC,B) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=TtrsB></algorithm>
        <output name=T> T </output>
      </compute>
      <var name=T operation=LE></var>
    </finish>
    <step>
      <var name=T><convert> T - 5 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=BCC_FCC></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2fcc> x(FCC,B) </input>
      <input once name=x2bcc> x(BCC,B) </input>
      <output> T </output>
      <output noname> x(FCC,B) </output>
      <output noname> x(BCC,B) </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output ChangeXY>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
        <output name=x2bcc> x(BCC_1,B) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=BCC_L_A2B></algorithm>
        <output name=T> T </output>
      </compute>
      <var name=T operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T + 5 </convert></var>
    </step>
    <compute>
      <PhaseEquilibrium class=algorithm id=BCC_A2B
        debug=0
        SaveSolution=0
        ThrowException=1>
        <phases>
          <phase class=phase IDREF=BCC></phase>
          <phase class=phase IDREF=A2B></phase>
        </phases>
        <state status=dependent> x(BCC,A) </state>
        <state status=unknown lower=0 upper=0.3333> x(BCC,B) </state>
        <state status=constraint> T </state>
        <state status=HardConstraint value=1> p </state>
      </PhaseEquilibrium> 
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2bcc> x(BCC,B) </input>
      <output> T </output>
      <output noname> x(BCC,B) </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output ChangeXY>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
        <output name=x2bcc> x(BCC_2,B) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=A2B_L_BCC></algorithm>
        <output name=T> T </output>
      </compute>
      <var name=T operation=GE></var>
    </finish>
    <step>
      <var name=T><convert> T + 5 </convert></var>
    </step>
    <compute>
      <PhaseEquilibrium class=algorithm id=A2B_BCC
        debug=0
        SaveSolution=0
        ThrowException=1>
        <phases>
          <phase class=phase IDREF=BCC></phase>
          <phase class=phase IDREF=A2B></phase>
        </phases>
        <state status=dependent> x(BCC,A) </state>
        <state status=unknown lower=0.33333 upper=1> x(BCC,B) </state>
        <state status=constraint> T </state>
        <state status=HardConstraint value=1> p </state>
      </PhaseEquilibrium> 
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2bcc> x(BCC,B) </input>
      <output> T </output>
      <output noname> x(BCC,B) </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output ChangeXY>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
        <output name=x2bcc1> x(BCC_1,B) </output>
        <output name=x2bcc2> x(BCC_2,B) </output>
        <output> T </output>
      </compute>
    </start>
    <finish>
      <var name=T value=500 operation=LE></var>
    </finish>
    <step>
      <var name=T><convert> T - 5 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=BCC_BCC></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=T> T </input>
      <input once name=x2bcc1> x(BCC_1,B) </input>
      <input once name=x2bcc2> x(BCC_2,B) </input>
      <output> T </output>
      <output noname> x(BCC_1,B) </output>
      <output noname> x(BCC_2,B) </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
        <output> T </output>
      </compute>
      <var name=x2 value=0.33333></var>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=A2B_L></algorithm>
        <input value=0.3333> x(L,B) </input>
        <output> T </output>
      </compute>
    </finish>
    <step>
      <var name=T><convert> T + 2000 </convert></var>
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
        <algorithm class=algorithm IDREF=BCC_L_A2B></algorithm>
        <output> T </output>
        <output name=x2> x(BCC,B) </output>
      </compute>
    </start>
    <finish>
      <var name=x2 value=0.33333 operation=GE></var>
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
      <compute> 
        <algorithm class=algorithm IDREF=A2B_L_BCC></algorithm>
        <output> T </output>
        <output name=x2> x(BCC,B) </output>
      </compute>
    </start>
    <finish>
      <var name=x2 value=0.33333 operation=LE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 - 1 </convert></var>
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
      <compute> 
        <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
        <output> T </output>
        <output name=x2> x(L,B) </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=L_BCC_FCC></algorithm>
        <output name=x2> x(FCC,B) </output>
      </compute>
      <var name=x2 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 + 1 </convert></var>
    </step>
    <compute> 
      <PassThrough class=algorithm id=p4></PassThrough>
      <input name=T> T </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
  <ComputeOutput class=output>
    <start>
      <compute> 
        <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
        <output> T </output>
        <output name=x2> x(BCC_1,B) </output>
      </compute>
    </start>
    <finish>
      <compute> 
        <algorithm class=algorithm IDREF=BCC_A2B_BCC></algorithm>
        <output name=x2> x(BCC_2,B) </output>
      </compute>
      <var name=x2 operation=GE></var>
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
</OutputFile>

<OutputFile ext=hmix format=gnuplot>
  <SeriesOutput Residuals=0> Hmix </SeriesOutput>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0.0></var>
    </start>
    <finish>
      <var name=x2 value=1.0 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2 + 0.01 </convert>
      </var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=L></algorithm>
      <input name=x2> x(B) </input>
      <input Value=1> p </input>
      <input Value=1773> T </input>
      <output> x(B) </output>
      <output> H_mix </output>
    </compute>
  </ComputeOutput>
</OutputFile>

<OutputFile ext=act format=gnuplot>
  <SeriesOutput Residuals=0> A </SeriesOutput>
  <ResidualOutput class=output InName=x2 OutName=a_B>
    <start>
      <var name=x2 value=0></var>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert>x2+0.01</convert></var>
    </step>
    act <var name=T value=1573></var>
  </ResidualOutput> 
</OutputFile>

<OutputFile ext=g600 format=gnuplot>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0></var>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=L></algorithm>
      <input once value=600> T </input>
      <input name=x2> x(2) </input>
      <output> x(2) </output>
      <output name=L> G </output>
    </compute>
    <compute>
      <algorithm IDREF=BCC></algorithm>
      <input once value=600> T </input>
      <input name=x2> x(2) </input>
      <output name=bcc> G </output>
    </compute>
    <compute>
      <algorithm IDREF=FCC></algorithm>
      <input once value=600> T </input>
      <input name=x2> x(2) </input>
      <output name=fcc> G </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <compute>
        <algorithm IDREF=A2B></algorithm>
        <input value=600> T </input>
        <output> G </output>
        <convert name=A2B> G/3 </convert>
      </compute>
      <var name=x2 value=0.3333></var>
    </start>
    <finish>
      <compute>
        <algorithm IDREF=A2B></algorithm>
        <input value=600> T </input>
        <output> G </output>
        <convert name=A2B> G/3 + 1000 </convert>
      </compute>
      <var name=A2B operation=GE></var>
    </finish>
    <step>
      <var name=A2B><convert> A2B+10000 </convert>
      </var>
    </step>
    <compute> 
      <algorithm class=algorithm IDREF=pass></algorithm>
      <input name=A2B> A2B </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output name=A2B> A2B </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

<OutputFile ext=g1100 format=gnuplot>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0></var>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=L></algorithm>
      <input once value=1100> T </input>
      <input name=x2> x(2) </input>
      <output> x(2) </output>
      <output name=L> G </output>
    </compute>
    <compute>
      <algorithm IDREF=BCC></algorithm>
      <input once value=1100> T </input>
      <input name=x2> x(2) </input>
      <output name=bcc> G </output>
    </compute>
    <compute>
      <algorithm IDREF=FCC></algorithm>
      <input once value=1100> T </input>
      <input name=x2> x(2) </input>
      <output name=fcc> G </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <compute>
        <algorithm IDREF=A2B></algorithm>
        <input value=1100> T </input>
        <output> G </output>
        <convert name=A2B> G/3 </convert>
      </compute>
      <var name=x2 value=0.3333></var>
    </start>
    <finish>
      <compute>
        <algorithm IDREF=A2B></algorithm>
        <input value=1100> T </input>
        <output> G </output>
        <convert name=A2B> G/3 + 1000 </convert>
      </compute>
      <var name=A2B operation=GE></var>
    </finish>
    <step>
      <var name=A2B><convert> A2B+10000 </convert></var>
    </step>
    <compute> 
      <algorithm class=algorithm IDREF=pass></algorithm>
      <input name=A2B> A2B </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output name=A2B> A2B </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

<OutputFile ext=g1500 format=gnuplot>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0></var>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.01 </convert></var>
    </step>
    <compute>
      <algorithm IDREF=L></algorithm>
      <input once value=1500> T </input>
      <input name=x2> x(2) </input>
      <output> x(2) </output>
      <output name=L> G </output>
    </compute>
    <compute>
      <algorithm IDREF=BCC></algorithm>
      <input once value=1500> T </input>
      <input name=x2> x(2) </input>
      <output name=bcc> G </output>
    </compute>
    <compute>
      <algorithm IDREF=FCC></algorithm>
      <input once value=1500> T </input>
      <input name=x2> x(2) </input>
      <output name=fcc> G </output>
    </compute>
  </ComputeOutput> 
  <ComputeOutput class=output>
    <start>
      <compute>
        <algorithm IDREF=A2B></algorithm>
        <input value=1500> T </input>
        <output> G </output>
        <convert name=A2B> G/3 </convert>
      </compute>
      <var name=x2 value=0.3333></var>
    </start>
    <finish>
      <compute>
        <algorithm IDREF=A2B></algorithm>
        <input value=1500> T </input>
        <output> G </output>
        <convert name=A2B> G/3 + 1000 </convert>
      </compute>
      <var name=A2B operation=GE></var>
    </finish>
    <step>
      <var name=A2B><convert> A2B+10000 </convert></var>
    </step>
    <compute> 
      <algorithm class=algorithm IDREF=pass></algorithm>
      <input name=A2B> A2B </input>
      <input name=x2> x2 </input>
      <output> x2 </output>
      <output name=A2B> A2B </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

