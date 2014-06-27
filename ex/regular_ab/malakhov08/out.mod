<OutputFile ext=pd format=gnuplot>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0.0></var>
    </start>
    <finish>
      <var name=x2 value=1.0 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.005 </convert></var>
    </step>
    <compute>
      <algorithm class=algorithm IDREF=A_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=x2> x(L,B) </input>
      <output> x(L,B) </output>
      <output noname> T </output>
    </compute>
  </ComputeOutput>
</OutputFile>

