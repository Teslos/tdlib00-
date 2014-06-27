<coef id=S> 1 </coef>
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

<PointPhase class=phase id=A>
  <species> A
    <Cp_zero class=func_Tp>
      (
				<coef computed id=H1> -S*Tm </coef> +
      	<coef IDREF=S> </coef> *T)
    </Cp_zero> 
  </species>
</PointPhase>

