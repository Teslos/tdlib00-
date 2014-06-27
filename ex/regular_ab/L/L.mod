<coef id=A> 0 </coef>
<elements>
  A    20
  B    50
</elements>
<SimpleSolution class=phase id=L>
  <components> A B </components>
  <Reference class=FuncTpx>
    <func_x> +x(A)* </func_x>
    <species>
      <Cp_zero class=func_Tp>
        (
        <coef> 0 </coef> +
        <coef> 0 </coef> *T)
      </Cp_zero> 
    </species>
    <func_x> +x(B)* </func_x>
    <species>
      <Cp_zero class=func_Tp>
        (
        <coef> 0 </coef> +
        <coef> 0 </coef> *T)
      </Cp_zero> 
    </species>
  </Reference>  
  <IdealMixing>
    +R*T*x(A)*log(x(A))
    +R*T*x(B)*log(x(B))
  </IdealMixing>
  <RedlichKister class=FuncTpx>
    +x(A)*x(B)*(
      <func_x> + </func_x> 
      <Cp_zero class=func_Tp>
        (
        <coef IDREF=A> </coef>+
        <coef> 0 </coef>*T)
      </Cp_zero> 
    )
  </RedlichKister> 
</SimpleSolution>


