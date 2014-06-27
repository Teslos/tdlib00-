<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=prop format=file>
  <ComputeOutput>
    <start>
      <var name=T value=950></var>
    </start>
    <finish>
      <var name=T value=850 operation=LE></var>
    </finish>
    <step>
      <var name=T><convert> T-25 </convert></var>
    </step>
    <compute>
      <PhaseEquilibrium class=algorithm id=s23_L_2>
        <phases>
          <phase class=phase IDREF=L></phase>
          <phase class=phase IDREF=s23></phase>
        </phases>
        <state status=dependent> x(L,Bi) </state>
        <state status=unknown> x(L,Se) </state>
        <state status=constraint> T </state>
        <state status=HardConstraint value=1> p </state>
      </PhaseEquilibrium>
      <input name=T> T </input>
      <output> T </output>
      <output> x(L,Se) </output>
    </compute>
  </ComputeOutput>
</OutputFile>
