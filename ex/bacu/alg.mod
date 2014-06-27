<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<PhaseProperty class=algorithm id=BaCu
   DependentMoleFraction=first>
   <phase class=phase IDREF=BaCu_s></phase>
</PhaseProperty> 
<PhaseProperty class=algorithm id=L>
  <phase class=phase IDREF=L></phase>
</PhaseProperty>
<PassThrough class=algorithm id=pass>
</PassThrough>

<PhaseEquilibrium class=algorithm id=BaCu13_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=BaCu13_s></phase>
  </phases>
  <state status=dependent> x(L,Ba) </state>
  <state status=constraint> x(L,Cu) </state>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=BaCu13_L_Cu
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=BaCu13_s></phase>
    <phase class=phase IDREF=Cu_s></phase>
  </phases>
  <state status=dependent> x(L,Ba) </state>
  <state status=unknown lower=0.5 upper=1> x(L,Cu) </state>
  <state status=unknown lower=900 upper=1100 value=1030> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=BaCu_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=BaCu_s></phase>
  </phases>
  <state status=dependent> x(L,Ba) </state>
  <state status=constraint> x(L,Cu) </state>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=BaCu_L_BaCu13
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=BaCu_s></phase>
    <phase class=phase IDREF=BaCu13_s></phase>
  </phases>
  <state status=dependent> x(L,Ba) </state>
  <state status=unknown lower=0.4 upper=0.6> x(L,Cu) </state>
  <state status=unknown lower=700 upper=900 value=820> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=Ba_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=Ba_s></phase>
  </phases>
  <state status=dependent> x(L,Ba) </state>
  <state status=constraint> x(L,Cu) </state>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=Ba_L_BaCu
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=Ba_s></phase>
    <phase class=phase IDREF=BaCu_s></phase>
  </phases>
  <state status=dependent> x(L,Ba) </state>
  <state status=unknown lower=0 upper=0.5> x(L,Cu) </state>
  <state status=unknown lower=600 upper=800 value=730> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=Cu_L
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=Cu_s></phase>
  </phases>
  <state status=dependent> x(L,Ba) </state>
  <state status=constraint> x(L,Cu) </state>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=TmBa
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=Ba_s></phase>
  </phases>
  <state status=dependent> x(L,Ba) </state>
  <state status=HardConstraint value=0> x(L,Cu) </state>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=TmCu
  debug=0
  SaveSolution=0
  ThrowException=1>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=Cu_s></phase>
  </phases>
  <state status=dependent> x(L,Ba) </state>
  <state status=HardConstraint value=1> x(L,Cu) </state>
  <state status=unknown> T </state>
  <state status=HardConstraint value=1> p </state>
</PhaseEquilibrium> 

