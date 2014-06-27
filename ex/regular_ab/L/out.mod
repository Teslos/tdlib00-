<PhaseProperty class=algorithm id=L
  DependentMoleFraction=first>
  <phase class=phase IDREF=L></phase>
</PhaseProperty> 

<OutputFile ext=G format=file>
  <ComputeOutput class=output>
    <start>
      <var name=x2 value=0></var>
    </start>
    <finish>
      <var name=x2 value=1 operation=GE></var>
    </finish>
    <step>
      <var name=x2><convert> x2+0.02 </convert>
      </var>
    </step>
    <compute>
      <algorithm IDREF=L></algorithm>
      <input value=1000> T </input>
      <input name=x2> x(B) </input>
      <output noname> x(B) </output>
      <output noname> G </output>
    </compute>
  </ComputeOutput> 
</OutputFile>

