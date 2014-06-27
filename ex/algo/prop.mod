<!--
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
-->

<OutputFile ext=prop format=file>
  <ComputeOutput>
    <compute>
      <PhaseEquilibrium class=algorithm id=prop>
        <phases>
          <PointPhase class=phase id=test>
            <species> K
              <Cp_const class=func_Tp>
                (
                <coef> 10000 </coef> +
                <coef> 5 </coef> *T+
                <coef> -1.5 </coef> *T*log(T))
              </Cp_const> 
            </species>
          </PointPhase>
        </phases>
        <state status=value value=6.39> S </state>
        <state status=unknown> T </state>
        <state status=HardConstraint value=1> p </state>
      </PhaseEquilibrium>
      <output> T </output>
      <output> S </output>
    </compute>
  </ComputeOutput>
</OutputFile>
