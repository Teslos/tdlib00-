<coef id=Tm> 2000 </coef>
<globals R=1 PETmax=2200 PETmin=10></globals>
<elements>
  A    20
  B    50
</elements>
<SimpleSolution class=phase id=L>
  <components> A B </components>
  <IdealMixing>
  </IdealMixing>
</SimpleSolution>

<PointPhase class=phase id=A1>
  <species> A
    <Cp_zero class=func_Tp>
      (
				<coef computed id=H1> -S1*Tm </coef> +
      	<coef id=S1> 0.25 </coef> *T)
    </Cp_zero> 
  </species>
</PointPhase>
<PointPhase class=phase id=A2>
  <species> A
    <Cp_zero class=func_Tp>
      (
				<coef computed id=H2> -S2*Tm </coef> +
      	<coef id=S2> 0.5 </coef> *T)
    </Cp_zero> 
  </species>
</PointPhase>
<PointPhase class=phase id=A3>
  <species> A
    <Cp_zero class=func_Tp>
      (
				<coef computed id=H3> -S3*Tm </coef> +
      	<coef id=S3> 1 </coef> *T)
    </Cp_zero> 
  </species>
</PointPhase>
<PointPhase class=phase id=A4>
  <species> A
    <Cp_zero class=func_Tp>
      (
				<coef computed id=H4> -S4*Tm </coef> +
      	<coef id=S4> 2 </coef> *T)
    </Cp_zero> 
  </species>
</PointPhase>
<PointPhase class=phase id=A5>
  <species> A
    <Cp_zero class=func_Tp>
      (
				<coef computed id=H5> -S5*Tm </coef> +
      	<coef id=S5> 4 </coef> *T)
    </Cp_zero> 
  </species>
</PointPhase>
<PointPhase class=phase id=A6>
  <species> A
    <Cp_zero class=func_Tp>
      (
				<coef computed id=H6> -S6*Tm </coef> +
      	<coef id=S6> 8 </coef> *T)
    </Cp_zero> 
  </species>
</PointPhase>

<PhaseEquilibrium class=algorithm id=A1_L>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=A1></phase>
  </phases>
  <state status=constraint> x(L,B) </state>
  <state status=unknown> T </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=A2_L>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=A2></phase>
  </phases>
  <state status=constraint> x(L,B) </state>
  <state status=unknown> T </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=A3_L>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=A3></phase>
  </phases>
  <state status=constraint> x(L,B) </state>
  <state status=unknown> T </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=A4_L>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=A4></phase>
  </phases>
  <state status=constraint> x(L,B) </state>
  <state status=unknown> T </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=A5_L>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=A5></phase>
  </phases>
  <state status=constraint> x(L,B) </state>
  <state status=unknown> T </state>
</PhaseEquilibrium> 
<PhaseEquilibrium class=algorithm id=A6_L>
  <phases>
    <phase class=phase IDREF=L></phase>
    <phase class=phase IDREF=A6></phase>
  </phases>
  <state status=constraint> x(L,B) </state>
  <state status=unknown> T </state>
</PhaseEquilibrium> 

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
      <algorithm class=algorithm IDREF=A1_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=x2> x(L,B) </input>
      <output> x(L,B) </output>
      <output name=0.25> T </output>
    </compute>
    <compute>
      <algorithm class=algorithm IDREF=A2_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=x2> x(L,B) </input>
      <output name=0.5> T </output>
    </compute>
    <compute>
      <algorithm class=algorithm IDREF=A3_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=x2> x(L,B) </input>
      <output name=1> T </output>
    </compute>
    <compute>
      <algorithm class=algorithm IDREF=A4_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=x2> x(L,B) </input>
      <output name=2> T </output>
    </compute>
    <compute>
      <algorithm class=algorithm IDREF=A5_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=x2> x(L,B) </input>
      <output name=4> T </output>
    </compute>
    <compute>
      <algorithm class=algorithm IDREF=A6_L></algorithm>
      <input once value=1> SaveSolution </input>
      <input name=x2> x(L,B) </input>
      <output name=8> T </output>
    </compute>
  </ComputeOutput>
</OutputFile>

