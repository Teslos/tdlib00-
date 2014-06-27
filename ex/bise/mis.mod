<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=prop format=file>
  <ComputeOutput>
    <start>
      <var name=T value=925></var>
    </start>
    <finish>
      <var name=T value=850 operation=LE></var>
    </finish>
    <step>
      <var name=T><convert> T-25 </convert></var>
    </step>
    <compute>
      <PhaseEquilibrium class=algorithm id=L1_L2>
        <phases>
          <phase class=phase IDREF=L></phase>
          <phase class=phase IDREF=L></phase>
        </phases>
        <state status=dependent> x(L_1,Bi) </state>
        <state status=unknown> x(L_1,Se) </state>
        <state status=dependent> x(L_2,Bi) </state>
        <state status=unknown> x(L_2,Se) </state>
        <state status=constraint> T </state>
        <state status=HardConstraint value=1> p </state> 
      </PhaseEquilibrium>
      <input name=T> T </input>
      <output> T </output>
      <output> x(L_1,Se) </output>
      <output> x(L_2,Se) </output>
    </compute>
  </ComputeOutput>
</OutputFile>
