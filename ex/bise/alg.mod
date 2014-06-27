<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<globals
  PETmin=300
  PETmax=1250>
</globals>
<PhaseEquilibrium class=algorithm id=L1_L2
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=L></phase>
  </phases>
  <state status=dependent> x(L_1,Bi) </state>
  <state status=unknown lower=0.6 upper=0.76 value=0.72> x(L_1,Se) </state>
  <state status=dependent> x(L_2,Bi) </state>
  <state status=unknown lower=0.89 upper=1 value=0.95> x(L_2,Se) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=Tm_s01
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s01></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=HardConstraint value=1> x(L,Se) </state>
  <state status=unknown lower=493 upper=494 value=493.9> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=Tm_s10
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s10></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=HardConstraint value=0> x(L,Se) </state>
  <state status=unknown lower=544 upper=545 value=544.55> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=Tm_s23
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s23></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=HardConstraint value=0.6> x(L,Se) </state>
  <state status=unknown upper=1250 value=990> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s10_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s10></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=constraint> x(L,Se) </state>
  <state status=unknown upper=1250> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s01_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s01></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=constraint> x(L,Se) </state>
  <state status=unknown upper=1250> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s10_L_s32
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s10></phase>
    <phase class=phase IDREF=s32></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=unknown lower=0 upper=0.1> x(L,Se) </state>
  <state status=unknown lower=500 upper=545 value=540> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s11_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s11></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=constraint> x(L,Se) </state>
  <state status=dependent> x(s11,Bi) </state>
  <state status=unknown lower=0.4 upper=0.6> x(s11,Se) </state>
  <state status=unknown upper=1250> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s11_L_1
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s11></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=unknown lower=0.15 upper=0.5> x(L,Se) </state>
  <state status=dependent> x(s11,Bi) </state>
  <state status=unknown lower=0.4 upper=0.6> x(s11,Se) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s11_L_s23
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s11></phase>
    <phase class=phase IDREF=s23></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=unknown lower=0.2 upper=0.6> x(L,Se) </state>
  <state status=dependent> x(s11,Bi) </state>
  <state status=unknown lower=0.4 upper=0.6> x(s11,Se) </state>
  <state status=unknown lower=800 upper=950 value=880> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s11_s23
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=s11></phase>
    <phase class=phase IDREF=s23></phase>
  </phases>
  <state status=dependent> x(s11,Bi) </state>
  <state status=unknown lower=0.4 upper=0.6> x(s11,Se) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s11_s32
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=s11></phase>
    <phase class=phase IDREF=s32></phase>
  </phases>
  <state status=dependent> x(s11,Bi) </state>
  <state status=unknown lower=0.4 upper=0.6> x(s11,Se) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s11_stoich
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=s11></phase>
  </phases>
  <state status=dependent> x(s11,Bi) </state>
  <state status=HardConstraint value=0.5> x(s11,Se) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s23_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s23></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=constraint> x(L,Se) </state>
  <state status=unknown upper=1250> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s23_L1_L2
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s23></phase>
  </phases>
  <state status=dependent> x(L_1,Bi) </state>
  <state status=unknown lower=0.66 upper=0.75 value=0.72> x(L_1,Se) </state>
  <state status=dependent> x(L_2,Bi) </state>
  <state status=unknown lower=0.92 upper=1 value=0.96> x(L_2,Se) </state>
  <state status=unknown lower=800 upper=1000 value=893> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s23_L_l
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s23></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=unknown lower=0.4 upper=0.6> x(L,Se) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s23_L_r
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s23></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=unknown lower=0.6 upper=1> x(L,Se) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s23_L_s01
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s23></phase>
    <phase class=phase IDREF=s01></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=unknown lower=0.96 upper=0.999> x(L,Se) </state>
  <state status=unknown lower=470 upper=493.84 value=489> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s32_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s32></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=constraint> x(L,Se) </state>
  <state status=unknown upper=1250> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s32_L_l
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s32></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=unknown lower=0 upper=0.2> x(L,Se) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=s32_L_s11
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=s32></phase>
    <phase class=phase IDREF=s11></phase>
  </phases>
  <state status=dependent> x(L,Bi) </state>
  <state status=unknown lower=0.15 upper=0.4> x(L,Se) </state>
  <state status=dependent> x(s11,Bi) </state>
  <state status=unknown lower=0.4 upper=0.6> x(s11,Se) </state>
  <state status=unknown lower=650 upper=850 value=744> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 

