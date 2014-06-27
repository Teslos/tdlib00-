<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PhaseProperty class=algorithm id=A2B
    DependentMoleFraction=first>
   <phase class=phase IDREF=A2B></phase>
</PhaseProperty> 
<PhaseProperty class=algorithm id=L
    DependentMoleFraction=first>
    <phase class=phase IDREF=L></phase>
</PhaseProperty>
<PhaseProperty class=algorithm id=Bfcc
    DependentMoleFraction=first>
    <phase class=phase IDREF=B_fcc></phase>
</PhaseProperty>
<PassThrough class=algorithm id=pass>
</PassThrough>

<PhaseEquilibrium class=algorithm id=BCC_L_l
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=BCC></phase>
  </phases>
  <state status=dependent> x(L,A) </state>
  <state status=unknown lower=0 upper=0.333333> x(L,B) </state>
  <state status=dependent> x(BCC,A) </state>
  <state status=unknown lower=0 upper=0.333333> x(BCC,B) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=BCC_L_r
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=BCC></phase>
  </phases>
  <state status=dependent> x(L,A) </state>
  <state status=unknown lower=0.33333 upper=1> x(L,B) </state>
  <state status=dependent> x(BCC,A) </state>
  <state status=unknown lower=0.33333 upper=1> x(BCC,B) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=BCC_BCC
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=BCC></phase>
    <phase class=phase IDREF=BCC></phase>
  </phases>
  <state status=dependent> x(BCC_1,A) </state>
  <state status=unknown lower=0 upper=0.2> x(BCC_1,B) </state>
  <state status=dependent> x(BCC_2,A) </state>
  <state status=unknown lower=0.7 upper=1> x(BCC_2,B) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=FCC_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=FCC></phase>
  </phases>
  <state status=dependent> x(L,A) </state>
  <state status=unknown lower=0.33333 upper=1> x(L,B) </state>
  <state status=dependent> x(FCC,A) </state>
  <state status=unknown lower=0.33333 upper=1> x(FCC,B) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=BCC_FCC
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=FCC></phase>
    <phase class=phase IDREF=BCC></phase>
  </phases>
  <state status=dependent> x(FCC,A) </state>
  <state status=unknown> x(FCC,B) </state>
  <state status=dependent> x(BCC,A) </state>
  <state status=unknown> x(BCC,B) </state>
  <state status=constraint> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=A2B_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=A2B></phase>
  </phases>
  <state status=dependent> x(L,A) </state>
  <state status=constraint> x(L,B) </state>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 

<PhaseEquilibrium class=algorithm id=BCC_L_A2B
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=BCC></phase>
    <phase class=phase IDREF=A2B></phase>
  </phases>
  <state status=dependent> x(L,A) </state>
  <state status=unknown lower=0 upper=0.3333 value=0.22> x(L,B) </state>
  <state status=dependent> x(BCC,A) </state>
  <state status=unknown lower=0 upper=0.15 value=0.056> x(BCC,B) </state>
  <state status=unknown lower=1000 upper=1300 value=1193> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=A2B_L_BCC
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=BCC></phase>
    <phase class=phase IDREF=A2B></phase>
  </phases>
  <state status=dependent> x(L,A) </state>
  <state status=unknown lower=0.33333 upper=1 value=0.52> x(L,B) </state>
  <state status=dependent> x(BCC,A) </state>
  <state status=unknown lower=0.33333 upper=1 value=0.80> x(BCC,B) </state>
  <state status=unknown lower=900 upper=1200 value=1049> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=L_BCC_FCC
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=BCC></phase>
    <phase class=phase IDREF=FCC></phase>
  </phases>
  <state status=dependent> x(L,A) </state>
  <state status=unknown lower=0.33333 upper=1 value=0.63> x(L,B) </state>
  <state status=dependent> x(BCC,A) </state>
  <state status=unknown lower=0.7 upper=1 value=0.84> x(BCC,B) </state>
  <state status=dependent> x(FCC,A) </state>
  <state status=unknown lower=0.7 upper=1 value=0.86> x(FCC,B) </state>
  <state status=unknown lower=1100 upper=1300 value=1203> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=BCC_A2B_BCC
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=BCC></phase>
    <phase class=phase IDREF=BCC></phase>
    <phase class=phase IDREF=A2B></phase>
  </phases>
  <state status=dependent> x(BCC_1,A) </state>
  <state status=unknown lower=0 upper=0.2 value=0.015> x(BCC_1,B) </state>
  <state status=dependent> x(BCC_2,A) </state>
  <state status=unknown lower=0.7 upper=1 value=0.76> x(BCC_2,B) </state>
  <state status=unknown lower=600 upper=800 value=726> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 

<PhaseEquilibrium class=algorithm id=TmA
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=A_bcc></phase>
    <phase class=phase IDREF=A_l></phase>
  </phases>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=TmB
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=B_fcc></phase>
    <phase class=phase IDREF=B_l></phase>
  </phases>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=TtrsB
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=B_bcc></phase>
    <phase class=phase IDREF=B_fcc></phase>
  </phases>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 


