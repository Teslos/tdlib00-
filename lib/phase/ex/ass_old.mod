<associated_solution class=phase id=L>
  <components> Hg Te </components>
  <Reference class=FuncTpx>
    <func_x> +x(Hg)* </func_x>
    <species id=Hg> Hg
    </species>
    <func_x> +x(Te)* </func_x>
    <species id=Te> Te
    </species>
  </Reference> 
  <internal_solution>
    <SimpleSolution class=phase>
      <components> Hg_l Te_l Hg0.5Te0.5_l </components>
      <Reference class=FuncTpx>
        <func_x> +x(Hg_l)* </func_x>
        <species id=Hg_l> Hg_l
        </species>
        <func_x> +x(Te_l)* </func_x>
        <species id=Te_l> Te_l
        </species>
        <func_x> +x(Hg0.5Te0.5_l)* </func_x>
        <species id=Hg0.5Te0.5_l> Hg0.5Te0.5_l
          <Cp_zero class=func_Tp>
            (
            <coef> -10000 </coef> +
            <coef> 0 </coef> *T)
          </Cp_zero> 
        </species>
      </Reference> 
      <IdealMixing class=FuncTpx>
        +R*T*x(Hg_l)*log(x(Hg_l))
        +R*T*x(Te_l)*log(x(Te_l))
        +R*T*x(Hg0.5Te0.5_l)*log(x(Hg0.5Te0.5_l))
      </IdealMixing> 
      <RedlichKister class=FuncTpx formalism=Muggianu>
        +x(Hg_l)*x(Hg0.5Te0.5_l)*(
          <func_x> + </func_x> 
          <Cp_zero class=func_Tp>
            (
            <coef> 0 </coef> +
            <coef> 0 </coef> *T)
          </Cp_zero> 
        )
      </RedlichKister> 
      <RedlichKister class=FuncTpx formalism=Muggianu>
        +x(Te_l)*x(Hg0.5Te0.5_l)*(
          <func_x> + </func_x> 
          <Cp_zero class=func_Tp>
            (
            <coef> 0 </coef> +
            <coef> 0 </coef> *T)
          </Cp_zero> 
        )
      </RedlichKister> 
    </SimpleSolution> 
  </internal_solution>
</associated_solution> 

